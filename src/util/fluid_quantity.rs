use crate::util::helper::{cubic_pulse, length};
use std::mem::swap;
use crate::solid::SolidBody;

pub struct FluidQuantity {
    pub src: Vec<f64>,
    pub dst: Vec<f64>,
    pub normal_x:  Vec<f64>,
    pub normal_y:  Vec<f64>,
    pub cell: Vec<u8>,
    pub body: Vec<u8>,
    pub mask: Vec<u8>,
    pub rows:      usize,
    pub columns:   usize,
    pub x_offset:  f64,
    pub y_offset:  f64,
    pub cell_size: f64,
}

impl FluidQuantity {
    pub fn new(rows: usize, columns: usize, x_offset: f64, y_offset: f64, cell_size: f64) -> FluidQuantity {
        FluidQuantity {
            src: vec![0.0; rows * columns],
            dst: vec![0.0; rows * columns],
            normal_x: vec![0.0; rows * columns],
            normal_y: vec![0.0; rows * columns],
            cell: vec![0u8; rows * columns],
            body: vec![0u8; rows * columns],
            mask: vec![0u8; rows * columns],
            rows,
            columns,
            x_offset,
            y_offset,
            cell_size,
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

    pub fn body_at(&self, row: usize, column: usize) -> u8 {
        self.body[row * self.columns + column]
    }

    pub fn body_at_mut(&mut self, row: usize, column: usize) -> &mut u8 {
        &mut self.body[row * self.columns + column]
    }

    pub fn cell_at(&self, row: usize, column: usize) -> u8 {
        self.cell[row * self.columns + column]
    }

    pub fn cell_at_mut(&mut self, row: usize, column: usize) -> &mut u8 {
        &mut self.cell[row * self.columns + column]
    }

    pub fn mask_at(&self, row: usize, column: usize) -> u8 {
        self.mask[row * self.columns + column]
    }

    pub fn mask_at_mut(&mut self, row: usize, column: usize) -> &mut u8 {
        &mut self.mask[row * self.columns + column]
    }

    pub fn normal_x_at(&self, row: usize, column: usize) -> f64 {
        self.normal_x[row * self.columns + column]
    }

    pub fn normal_x_at_mut(&mut self, row: usize, column: usize) -> &mut f64 {
        &mut self.normal_x[row * self.columns + column]
    }

    pub fn normal_y_at(&self, row: usize, column: usize) -> f64 {
        self.normal_y[row * self.columns + column]
    }

    pub fn normal_y_at_mut(&mut self, row: usize, column: usize) -> &mut f64 {
        &mut self.normal_y[row * self.columns + column]
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

    pub fn fill_solid_fields(&mut self, bodies: &Vec<SolidBody>) {
        if bodies.is_empty() {
            return;
        }

        for row in 0..self.rows {
            for column in 0..self.columns {
                let x = (column as f64 + self.x_offset) * self.cell_size;
                let y = (row as f64 + self.y_offset) * self.cell_size;

                *self.body_at_mut(row, column) = 0;
                let mut d = bodies[0].distance(x, y);
                for index in 1..bodies.len() {
                    let id = bodies[index].distance(x, y);

                    if id < d {
                        // Note here that if there are more than 256 objects there will be errors
                        *self.body_at_mut(row, column) = index as u8;
                        d = id;
                    }
                }

                if d < 0.0 {
                    *self.cell_at_mut(row, column) = 1;
                } else {
                    *self.cell_at_mut(row, column) = 0;
                }

                // This needs to be refactored
                let output = bodies[self.body_at(row, column) as usize].distance_normal(x, y);
                *self.normal_x_at_mut(row, column) = output.0;
                *self.normal_y_at_mut(row, column) = output.1;
            }
        }
    }

    pub fn fill_solid_mask(&mut self) {
        for row in 0..(self.rows - 1) {
            for column in 0..(self.columns - 1) {
                if self.cell_at(row, column) == 1 {
                    let nx = self.normal_x_at(row, column);
                    let ny = self.normal_y_at(row, column);

                    *self.mask_at_mut(row, column) = 0;
                    if nx != 0.0 && self.cell_at(row, column + nx.signum() as usize) != 0 {
                        *self.mask_at_mut(row, column) |= 1;
                    }

                    if ny != 0.0 && self.cell_at(row + ny.signum() as usize, column) != 0 {
                        *self.mask_at_mut(row, column) |= 2;
                    }
                }
            }
        }
    }

    pub fn extrapolate_normal(&mut self, idx: usize) -> f64 {
        let nx = self.normal_x[idx];
        let ny = self.normal_y[idx];

        let src_x = self.src[idx +  nx.signum() as usize];
        let src_y = self.src[idx + ny.signum() as usize * self.columns];

        (nx.abs() * src_x + ny.abs() * src_y) / (nx.abs() + ny.abs())
    }

    pub fn free_neighbour(&mut self, idx: usize, border: &mut Vec<usize>, mask: u8) {
        self.mask[idx] &= !mask;
        if self.cell[idx] != 0 && self.mask[idx] == 0 {
            border.push(idx);
        }
    }

    pub fn extrapolate(&mut self) {
        self.fill_solid_mask();

        let mut border = vec![0 as usize; 0];

        for row in 0..(self.rows - 1) {
            for column in 0..(self.columns - 1) {
                let idx = column + row * self.columns;

                if self.cell[idx] != 0 && self.mask[idx] == 0 {
                    border.push(idx);
                }
            }
        }

        while !border.is_empty() {
            let idx = border.pop().unwrap();

            self.src[idx] = self.extrapolate_normal(idx);

            if self.normal_x[idx - 1] > 0.0 {
                self.free_neighbour(idx - 1, &mut border, 1);
            }
            if self.normal_x[idx + 1] < 0.0 {
                self.free_neighbour(idx + 1, &mut border, 1);
            }
            if self.normal_y[idx - self.columns] > 0.0 {
                self.free_neighbour(idx - self.columns, &mut border, 2);
            }
            if self.normal_y[idx + self.columns] < 0.0 {
                self.free_neighbour(idx + self.columns, &mut border, 2);
            }
        }
    }
}