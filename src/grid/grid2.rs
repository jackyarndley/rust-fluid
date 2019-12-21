pub trait Grid2<T> {
    fn test() {
        println!("This is a test!");
    }

    fn at(self, i: usize, j: usize) -> T {
        self.data[i * self.columns + j]
    }

    fn at_mut(&mut self, i: usize, j: usize) -> &mut T {
        &mut self.data[i * self.columns + j]
    }
}