use crate::util::{cubic_pulse, length};

pub struct Field {
    pub field: Vec<f64>,
    //pub field_dst: Vec<f64>,
    pub rows:      usize,
    pub columns:   usize,
    pub x_offset:  f64,
    pub y_offset:  f64
}

impl Field {
    pub fn new(rows: usize, columns: usize, x_offset: f64, y_offset: f64) -> Field {
        Field {
            field: vec![0.0; rows * columns],
            //field_dst: vec![0.0; rows * columns],
            rows,
            columns,
            x_offset,
            y_offset
        }
    }

    pub fn at(&self, row: usize, column: usize) -> f64 {
        self.field[row * self.columns + column]
    }

    pub fn at_mut(&mut self, row: usize, column: usize) -> &mut f64 {
        &mut self.field[row * self.columns + column]
    }

    pub fn add_inflow(&mut self, ix0: usize, iy0: usize, ix1: usize, iy1: usize, value: f64) {
        for y in iy0..iy1 {
            for x in ix0..ix1 {
                let l = length(
                    (2.0 * (x as f64) - (ix0 + ix1) as f64) / ((ix1 - ix0) as f64),
                    (2.0 * (y as f64) - (iy0 + iy1) as f64) / ((iy1 - iy0) as f64)
                );

                let vi = cubic_pulse(l) * value;

                if self.field[x + y * self.columns].abs() < vi.abs() {
                    self.field[x + y * self.columns] = vi;
                }
            }
        }
    }
}