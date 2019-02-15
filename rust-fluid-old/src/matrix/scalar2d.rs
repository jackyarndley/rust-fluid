pub struct Scalar2D<T> {
    pub array: Vec<T>,
    pub rows: usize,
    pub columns: usize,
}

impl<T> Scalar2D<T>{
    pub fn new(width: usize, height: usize) -> Scalar2D<T> {
        Scalar2D {
            array: vec![T; rows * columns],
            rows,
            columns
        }
    }

    pub fn at(&self, row: usize, column: usize) -> T {
        &self.array[self.columns * row + column]
    }

    pub fn at_mut(&mut self, row: usize, column: usize) -> &mut T {
        &mut self.array[self.columns * row + column]
    }
}