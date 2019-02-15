pub struct Vector2D<T> {
    pub array: Vec<T>,
    pub width: usize,
    pub height: usize,
    pub depth: usize,
}

impl<T> Vector2D<T>{
    pub fn new(width: usize, height: usize, depth: usize) -> Vector2D<T> {
        Vector2D {
            array: vec![T; width * height * depth],
            width,
            height,
            depth
        }
    }

    pub fn at(&self, x: usize, y: usize, z: usize) -> &T {
        &self.array[self.depth * (self.height * y + x) + z]
    }

    pub fn at_mut(&mut self, x: usize, y: usize, z: usize) -> &mut T {
        &mut self.array[self.depth * (self.height * y + x) + z]
    }
}