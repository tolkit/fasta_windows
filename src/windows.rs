
// clone of https://stackoverflow.com/questions/51257304/creating-a-sliding-window-iterator-of-slices-of-chars-from-a-string

// step_by should be the same as win_size in the generation of windows
// step_by should == 1 when computing kmer frequencies.

// TODO: is this efficient?

pub mod windows {

    pub fn char_windows<'a>(src: &'a str, win_size: usize, step_by: usize) -> impl Iterator<Item = &'a str> {
        src.char_indices()
            .step_by(step_by) // EDIT: skip the intervening windows, otherwise windows overlap
            .flat_map(move |(from, _)| {
                src[from ..].char_indices()
                    .skip(win_size - 1)
                    .next()
                    .map(|(to, c)| {
                        &src[from .. from + to + c.len_utf8()]
                    })
        })
    }

}