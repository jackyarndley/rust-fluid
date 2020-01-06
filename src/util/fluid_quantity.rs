use crate::util::helper::{cubic_pulse, length};
use std::mem::swap;

pub struct FluidQuantity {
    pub src: Vec<f64>,
    pub dst: Vec<f64>,
    pub rows:      usize,
    pub columns:   usize,
    pub x_offset:  f64,
    pub y_offset:  f64,
}

impl FluidQuantity {
    pub fn new(rows: usize, columns: usize, x_offset: f64, y_offset: f64) -> FluidQuantity {
        FluidQuantity {
            src: vec![0.0; rows * columns],
            dst: vec![0.0; rows * columns],
            rows,
            columns,
            x_offset,
            y_offset,
        }
    }

    pub fn at(&self, row: usize, column: usize) -> f64 {
        self.src[row * self.columns + column]
    }

    pub fn at_mut(&mut self, row: usize, column: usize) -> &mut f64 {
        &mut self.src[row * self.columns + column]
    }

    pub fn dst_at(&self, row: usize, column: usize) -> f64 {
        self.dst[row * self.columns + column]
    }

    pub fn dst_at_mut(&mut self, row: usize, column: usize) -> &mut f64 {
        &mut self.dst[row * self.columns + column]
    }

    pub fn swap_buffers(&mut self) {
        swap(&mut self.src, &mut self.dst);
    }

    pub fn add_inflow(&mut self, ix0: usize, iy0: usize, ix1: usize, iy1: usize, value: f64) {
        for y in iy0..iy1 {
            for x in ix0..ix1 {
                let l = length(
                    (2.0 * (x as f64) - (ix0 + ix1) as f64) / ((ix1 - ix0) as f64),
                    (2.0 * (y as f64) - (iy0 + iy1) as f64) / ((iy1 - iy0) as f64)
                );

                let vi = cubic_pulse(l) * value;

                if self.src[x + y * self.columns].abs() < vi.abs() {
                    self.src[x + y * self.columns] = vi;
                }
            }
        }
    }
}