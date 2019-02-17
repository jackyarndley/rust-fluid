use crate::interpolation::min;

pub struct Field {
    pub field:    Vec<f64>,
    pub rows:     usize,
    pub columns:  usize,
    pub offset_x: f64,
    pub offset_y: f64
}

impl Field {
    pub fn new(rows: usize, columns: usize, offset_x: f64, offset_y: f64) -> Field {
        Field {
            field: vec![0.0; rows * columns],
            rows,
            columns,
            offset_x,
            offset_y
        }
    }

    pub fn at(&self, row: usize, column: usize) -> f64 {
        self.field[row * self.columns + column]
    }

    pub fn at_mut(&mut self, row: usize, column: usize) -> &mut f64 {
        &mut self.field[row * self.columns + column]
    }

    // Returns length of vector
    fn length(x: f64, y: f64) -> f64 {
        (x * x + y * y).sqrt()
    }

    // Cubic pulse function returns value in range 0 - 1
    fn cubic_pulse(mut x: f64) -> f64 {
        x = min(x.abs(), 1.0);
        1.0 - x * x * (3.0 - 2.0 * x)
    }

    pub fn add_inflow(&mut self, ix0: usize, iy0: usize, ix1: usize, iy1: usize, value: f64) {
        for y in iy0..iy1 {
            for x in ix0..iy1 {
                let l = Field::length(
                    (2.0 * (x as f64) - (ix0 + ix1) as f64) / ((ix1 - ix0) as f64),
                    (2.0 * (y as f64) - (iy0 + iy1) as f64) / ((iy1 - iy0) as f64)
                );

                let vi = Field::cubic_pulse(l) * value;

                if self.field[x + y * self.columns].abs() < vi.abs() {
                    self.field[x + y * self.columns] = vi;
                }
            }
        }
    }
}