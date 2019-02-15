use std::mem::swap;

pub struct FluidQuantity {
    pub src: Vec<f64>,
    pub dst: Vec<f64>,
    pub width: usize,
    pub height: usize,
    pub offset_x: f64,
    pub offset_y: f64,
    pub cell_size: f64
}

impl FluidQuantity {
    pub fn new(width: usize, height: usize, offset_x: f64, offset_y: f64, cell_size: f64) -> FluidQuantity {
        FluidQuantity {
            src: vec![0f64; width * height],
            dst: vec![0f64; width * height],
            width,
            height,
            offset_x,
            offset_y,
            cell_size
        }
    }

    pub fn flip(&mut self) {
        swap(&mut self.src, &mut self.dst);
    }

    pub fn at(&self, x: usize, y: usize) -> f64 {
        self.src[self.width * y + x]
    }

    pub fn at_mut(&mut self, x: usize, y: usize) -> &mut f64 {
        &mut self.src[self.width * y + x]
    }

    pub fn grid_lerp(&self, x: f64, y: f64) -> f64 {

    }


}