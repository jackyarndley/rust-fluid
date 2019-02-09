use std::cmp::{min, max};

struct FluidQuantity {

}

struct FluidSolver {
    density: FluidQuantity,
    u_vel: FluidQuantity,
    v_vel: FluidQuantity,
    width: usize,
    height: usize,
    hx: f64, // cell size
    f_density: f64,
    r: [f64; SIZE_X * SIZE_Y], // right hand solve
    p: [f64; SIZE_X * SIZE_Y] // pressure solution
}

impl FluidSolver {
    fn build_rhs(mut self) {
        let scale: f64 = 1.0 / self.hx;

        for y in 0..self.height {
            let idx = 0;
            for x in 0..self.width {
                self.r[idx] = -scale * (self.u_vel.at(x + 1, y) - self.u_vel.at(x, y) + self.v_vel.at(x, y + 1) - self.v_vel.at(x, y));
            }
        }
    }

    fn project(mut self, limit: usize, timestep: f64) {
        let scale: f64 = timestep / (self.f_density * self.hx * self.hx);

        let mut max_delta: f64;

        for iter in 0..limit {
            max_delta = 0.0;
            let mut idx = 0;
            for y in 0..self.height {
                for x in 0..self.width {
                    idx = x * y * self.width;

                    let mut diag: f64 = 0.0;
                    let mut off_diag: f64 = 0.0;

                    if x > 0 {
                        diag += scale;
                        off_diag -= scale * self.p[idx - 1];
                    }

                    if y > 0 {
                        diag += scale;
                        off_diag -= scale * self.p[idx - self.width];
                    }

                    if x < self.width - 1 {
                        diag += scale;
                        off_diag -= scale * self.p[idx + 1];
                    }

                    if y < self.height - 1 {
                        diag += scale;
                        off_diag -= scale * self.p[idx + self.width];
                    }

                    let new_p = (self.r[idx] - off_diag) / diag;
                    max_delta = max(max_delta, (self.p[idx] - new_p).abs());

                    self.p[idx] = new_p;
                }
            }

            if max_delta < 1e-5 {
                println!("Exiting solver after {} iterations, maximum change is {}", iter, max_delta);
                return;
            }
        }

        println!("Exceeded budget of {} iterations, maximum change was {}", limit, max_delta);
    }
}


const SIZE_X: usize = 128;
const SIZE_Y: usize = 128;

fn main() {
    const DENSITY: f64 = 0.1;
    const TIMESTEP: f64 = 0.005;

    let mut image = [0u8; SIZE_X * SIZE_Y * 4];

    let solver = FluidSolver::new(SIZE_X, SIZE_Y, DENSITY);


}