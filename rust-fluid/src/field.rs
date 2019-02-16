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
}