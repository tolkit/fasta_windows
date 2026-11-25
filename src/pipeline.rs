//! Streaming FASTA → parallel window processing → ordered output.
//!
//! Three stages run concurrently:
//!   1. a reader thread parses the input into batches of whole records,
//!   2. rayon workers process the windows of the current batch, and
//!   3. a writer thread writes the previous batch's output.
//!
//! Each window's output is formatted straight into byte buffers inside the
//! worker, and results come back in input order, so nothing needs sorting and
//! memory is bounded by the batch size rather than the genome size.

use std::fs::File;
use std::io::Write;
use std::path::Path;
use std::sync::mpsc::{sync_channel, SyncSender};

use anyhow::{Context, Result};
use indicatif::{ProgressBar, ProgressStyle};
use needletail::errors::ParseErrorKind;
use needletail::parse_fastx_reader;
use rayon::prelude::*;

/// Bases of sequence collected before a batch is handed to the workers.
/// Records are never split, so a record longer than this is a batch on its own.
pub const BATCH_BASES: usize = 64 << 20;

/// Approximate amount of sequence per parallel task (rounded to whole windows).
pub const TASK_BASES: usize = 256 << 10;

pub struct Record {
    pub id: String,
    pub desc: Option<String>,
    pub seq: Vec<u8>,
}

/// One window of one record, as passed to the formatting callback.
pub struct Window<'a> {
    pub record: &'a Record,
    pub start: usize,
    pub end: usize,
    pub seq: &'a [u8],
}

#[inline]
fn is_header_ws(b: u8) -> bool {
    matches!(b, b' ' | b'\t' | b'\n' | b'\r' | 0x0b | 0x0c)
}

/// Split a FASTA header (without the leading '>') into ID and description, the
/// way bio::io::fasta did: trailing whitespace is trimmed, the ID runs up to the
/// first whitespace character, and the description is everything after it.
pub fn split_header(header: &[u8]) -> (String, Option<String>) {
    let end = header
        .iter()
        .rposition(|&b| !is_header_ws(b))
        .map_or(0, |i| i + 1);
    let header = &header[..end];
    match header.iter().position(|&b| is_header_ws(b)) {
        Some(i) => (
            String::from_utf8_lossy(&header[..i]).into_owned(),
            Some(String::from_utf8_lossy(&header[i + 1..]).into_owned()),
        ),
        None => (String::from_utf8_lossy(header).into_owned(), None),
    }
}

/// Stream `path` through `format` with the default batch and task sizes.
/// See [`run_with`].
pub fn run<W, F>(path: &Path, window_size: usize, writers: Vec<W>, format: F) -> Result<Vec<W>>
where
    W: Write + Send,
    F: Fn(&Window, &mut [Vec<u8>]) + Sync,
{
    run_with(path, window_size, BATCH_BASES, TASK_BASES, writers, format)
}

/// Stream every window of every record in `path` through `format`, in parallel,
/// and write the results to `writers` in input order.
///
/// `format` is called once per window and appends that window's text for
/// output `i` to `bufs[i]` (`bufs.len() == writers.len()`). Windows are
/// consecutive non-overlapping chunks of `window_size` bases; the last window
/// of a record may be shorter. The writers are returned unflushed.
pub fn run_with<W, F>(
    path: &Path,
    window_size: usize,
    batch_bases: usize,
    task_bases: usize,
    writers: Vec<W>,
    format: F,
) -> Result<Vec<W>>
where
    W: Write + Send,
    F: Fn(&Window, &mut [Vec<u8>]) + Sync,
{
    anyhow::ensure!(window_size > 0, "window size must be greater than zero");
    let file = File::open(path).with_context(|| format!("could not open {}", path.display()))?;
    let total_bytes = file.metadata().map(|m| m.len()).unwrap_or(0);
    let progress = ProgressBar::new(total_bytes);
    progress.set_style(
        ProgressStyle::with_template(
            "[+]\tReading fasta: {bar:40.cyan/blue} {bytes:>10}/{total_bytes:10}",
        )?
        .progress_chars(">>-"),
    );
    let n_out = writers.len();
    let windows_per_task = (task_bases / window_size).max(1);

    let writers = std::thread::scope(|scope| -> Result<Vec<W>> {
        // Rendezvous channels: the reader parses one batch ahead and the writer
        // writes one batch behind, no further, so at most two batches of
        // sequence and two of output are in memory at once.
        let (batch_tx, batch_rx) = sync_channel::<Result<Vec<Record>>>(0);
        let (out_tx, out_rx) = sync_channel::<Vec<Vec<Vec<u8>>>>(0);

        let reader_progress = progress.clone();
        scope.spawn(move || read_batches(file, reader_progress, batch_bases, batch_tx));

        let writer = scope.spawn(move || -> Result<Vec<W>> {
            let mut writers = writers;
            for chunks in out_rx {
                for bufs in chunks {
                    for (w, buf) in writers.iter_mut().zip(&bufs) {
                        w.write_all(buf)?;
                    }
                }
            }
            Ok(writers)
        });

        let mut result = Ok(());
        for batch in batch_rx {
            let batch = match batch {
                Ok(b) => b,
                Err(e) => {
                    result = Err(e);
                    break;
                }
            };
            let out = process_batch(&batch, window_size, windows_per_task, n_out, &format);
            if out_tx.send(out).is_err() {
                break; // writer failed; its error is reported by join below
            }
        }
        drop(out_tx);
        let writers = writer.join().expect("writer thread panicked")?;
        result.with_context(|| format!("failed to read {}", path.display()))?;
        Ok(writers)
    })?;
    progress.finish();
    Ok(writers)
}

fn read_batches(
    file: File,
    progress: ProgressBar,
    batch_bases: usize,
    tx: SyncSender<Result<Vec<Record>>>,
) {
    let mut reader = match parse_fastx_reader(progress.wrap_read(file)) {
        Ok(r) => r,
        Err(e) if e.kind == ParseErrorKind::EmptyFile => return,
        Err(e) => {
            let _ = tx.send(Err(e.into()));
            return;
        }
    };
    let mut batch = Vec::new();
    let mut bases = 0;
    while let Some(rec) = reader.next() {
        let rec = match rec {
            Ok(r) => r,
            Err(e) => {
                let _ = tx.send(Err(e.into()));
                return;
            }
        };
        let (id, desc) = split_header(rec.id());
        let seq = rec.seq().into_owned();
        bases += seq.len();
        batch.push(Record { id, desc, seq });
        if bases >= batch_bases {
            if tx.send(Ok(std::mem::take(&mut batch))).is_err() {
                return;
            }
            bases = 0;
        }
    }
    if !batch.is_empty() {
        let _ = tx.send(Ok(batch));
    }
}

/// Format every window in `batch`. Windows from all records are numbered
/// consecutively and cut into tasks of `windows_per_task`, so a large
/// chromosome is spread over many tasks and many small scaffolds share one.
/// The indexed parallel collect keeps tasks, and so windows, in input order.
fn process_batch<F>(
    batch: &[Record],
    window_size: usize,
    windows_per_task: usize,
    n_out: usize,
    format: &F,
) -> Vec<Vec<Vec<u8>>>
where
    F: Fn(&Window, &mut [Vec<u8>]) + Sync,
{
    // offsets[r] = number of windows in records before r
    let mut offsets = Vec::with_capacity(batch.len() + 1);
    offsets.push(0usize);
    for r in batch {
        offsets.push(offsets[offsets.len() - 1] + r.seq.len().div_ceil(window_size));
    }
    let total = offsets[batch.len()];

    (0..total.div_ceil(windows_per_task))
        .into_par_iter()
        .map(|task| {
            let first = task * windows_per_task;
            let last = (first + windows_per_task).min(total);
            let mut bufs = vec![Vec::new(); n_out];
            // last record whose first window is <= `first` (skips empty records)
            let mut r = offsets.partition_point(|&o| o <= first) - 1;
            for g in first..last {
                while offsets[r + 1] <= g {
                    r += 1;
                }
                let record = &batch[r];
                let start = (g - offsets[r]) * window_size;
                let end = (start + window_size).min(record.seq.len());
                let window = Window {
                    record,
                    start,
                    end,
                    seq: &record.seq[start..end],
                };
                format(&window, &mut bufs);
                if g == first {
                    // size the buffers from the first window, with headroom,
                    // rather than letting them double as they grow
                    let remaining = last - first - 1;
                    for buf in bufs.iter_mut() {
                        let per_window = buf.len() + buf.len() / 4;
                        buf.reserve_exact(per_window * remaining);
                    }
                }
            }
            bufs
        })
        .collect()
}

/// Append the decimal representation of `n`.
#[inline]
pub fn push_uint(buf: &mut Vec<u8>, mut n: u64) {
    let mut tmp = [0u8; 20];
    let mut i = tmp.len();
    loop {
        i -= 1;
        tmp[i] = b'0' + (n % 10) as u8;
        n /= 10;
        if n == 0 {
            break;
        }
    }
    buf.extend_from_slice(&tmp[i..]);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_split_header() {
        assert_eq!(split_header(b"chr1"), ("chr1".into(), None));
        assert_eq!(split_header(b"chr1 "), ("chr1".into(), None));
        assert_eq!(
            split_header(b"chr1 some description\r"),
            ("chr1".into(), Some("some description".into()))
        );
        assert_eq!(
            split_header(b"chr1\tsome  description"),
            ("chr1".into(), Some("some  description".into()))
        );
        // like bio: only the first whitespace character is consumed
        assert_eq!(split_header(b"chr1  x"), ("chr1".into(), Some(" x".into())));
        assert_eq!(split_header(b""), ("".into(), None));
    }

    #[test]
    fn test_push_uint() {
        for n in [0u64, 7, 10, 99, 100, 123456789, u64::MAX] {
            let mut b = Vec::new();
            push_uint(&mut b, n);
            assert_eq!(b, n.to_string().into_bytes());
        }
    }

    /// Batching and task splitting must not change which windows are produced
    /// or their order, including empty records, records shorter than a window,
    /// and tasks that span several records.
    #[test]
    fn test_window_order_across_batches_and_tasks() {
        let fasta =
            b">a desc a\nACGTACGTAC\nGT\n>empty\n\n>b\nAC\n>c\nACGTACGTACGTACGTACGTA\n>d x\nA\n";
        let dir = std::env::temp_dir().join(format!("fw_pipeline_test_{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("t.fa");
        std::fs::write(&path, fasta).unwrap();

        let fmt = |w: &Window, bufs: &mut [Vec<u8>]| {
            let line = format!(
                "{}\t{}\t{}\t{}\t{}\n",
                w.record.id,
                w.record.desc.as_deref().unwrap_or("-"),
                w.start,
                w.end,
                String::from_utf8_lossy(w.seq)
            );
            bufs[0].extend_from_slice(line.as_bytes());
            bufs[1].extend_from_slice(w.record.id.as_bytes());
        };
        let reference = String::from_utf8(
            run_with(
                &path,
                5,
                usize::MAX,
                usize::MAX,
                vec![Vec::new(), Vec::new()],
                fmt,
            )
            .unwrap()
            .remove(0),
        )
        .unwrap();
        assert_eq!(
            reference,
            "a\tdesc a\t0\t5\tACGTA\na\tdesc a\t5\t10\tCGTAC\na\tdesc a\t10\t12\tGT\n\
             b\t-\t0\t2\tAC\n\
             c\t-\t0\t5\tACGTA\nc\t-\t5\t10\tCGTAC\nc\t-\t10\t15\tGTACG\nc\t-\t15\t20\tTACGT\nc\t-\t20\t21\tA\n\
             d\tx\t0\t1\tA\n"
        );
        for batch_bases in [1, 3, 12, 15, 40] {
            for task_bases in [1, 5, 9, 10, 20, 1000] {
                let out = run_with(
                    &path,
                    5,
                    batch_bases,
                    task_bases,
                    vec![Vec::new(), Vec::new()],
                    fmt,
                )
                .unwrap();
                assert_eq!(String::from_utf8(out[0].clone()).unwrap(), reference);
                assert_eq!(out[1], b"aaabcccccd");
            }
        }

        std::fs::write(&path, b"").unwrap();
        let out = run_with(&path, 5, 10, 10, vec![Vec::<u8>::new()], fmt).unwrap();
        assert!(out[0].is_empty());
        std::fs::remove_dir_all(&dir).unwrap();
    }
}
