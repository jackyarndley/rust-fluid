use std::cmp::{min, max};
use std::mem::swap;
use std::io::{Write, stdout};
use image::save_buffer;

struct FluidQuantity {
    src: Vec<f64>,
    dst: Vec<f64>,
    width: usize,
    height: usize,
    offset_x: f64,
    offset_y: f64,
    cell_h: f64
}

impl FluidQuantity {
    pub fn new(width: usize, height: usize, offset_x: f64, offset_y: f64, cell_h: f64) -> FluidQuantity {
        FluidQuantity {
            src: vec!(0.0; width * height),
            dst: vec!(0.0; width * height),
            width,
            height,
            offset_x,
            offset_y,
            cell_h
        }
    }

    pub fn flip(&mut self) {
        swap(&mut self.src, &mut self.dst)
    }

    pub fn at(&self, x: usize, y: usize) -> f64 {
        self.src[x + y * self.width]
    }

    pub fn lerp(&self, x: f64, y: f64) -> f64 {
        0.0
    }

    pub fn advect(&self) {

    }

    pub fn add_inflow(&self) {

    }
}

struct FluidSolver {
    density: FluidQuantity,
    u_vel: FluidQuantity,
    v_vel: FluidQuantity,
    width: usize,
    height: usize,
    hx: f64, // cell size
    f_density: f64,
    r: Box<[f64; SIZE_X * SIZE_Y]>, // right hand solve
    p: Box<[f64; SIZE_X * SIZE_Y]> // pressure solution
}

impl FluidSolver {
    fn build_rhs(&mut self) {
        let scale: f64 = 1.0 / self.hx;

        for y in 0..self.height {
            let idx = 0;
            for x in 0..self.width {
                self.r[idx] = -scale * (self.u_vel.at(x + 1, y) - self.u_vel.at(x, y) + self.v_vel.at(x, y + 1) - self.v_vel.at(x, y));
            }
        }
    }

    fn project(&mut self, limit: usize, timestep: f64) {
        let scale: f64 = timestep / (self.f_density * self.hx * self.hx);

        let mut max_delta: f64 = 0.0;

        for iter in 0..limit {
            max_delta = 0.0;
            let mut idx = 0;
            for y in 0..self.height {
                for x in 0..self.width {
                    idx = x + y * self.width;

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

                    let temp = (self.p[idx] - new_p).abs();

                    if max_delta < temp {
                        max_delta = temp;
                    }

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

    fn apply_pressure(&self, timestep: f64) {
        let scale: f64 = timestep / (self.f_density * self.hx);

        let mut idx = 0;
        for y in 0..self.height {
            for x in 0..self.width {
                //self.u.at


                idx += 1;
            }
        }

        for y in 0..self.height {
            //
        }

        for x in 0..self.width {
            //
        }
    }

    pub fn new(width: usize, height: usize, density: f64) -> FluidSolver {
        let hx = 1.0 / min(width, height) as f64;
        FluidSolver {
            density: FluidQuantity::new(width, height, 0.5, 0.5, hx),
            u_vel: FluidQuantity::new(width + 1, height, 0.0, 0.5, hx),
            v_vel: FluidQuantity::new(width, height + 1, 0.5, 0.0, hx),
            width,
            height,
            hx,
            f_density: density,
            r: Box::new([0.0; SIZE_X * SIZE_Y]),
            p: Box::new([0.0; SIZE_X * SIZE_Y])
        }
    }

    pub fn update(&mut self, timestep: f64) {
        self.build_rhs();
        self.project(600, timestep);
        self.apply_pressure(timestep);

        self.density.advect();
        self.u_vel.advect();
        self.v_vel.advect();

        self.density.flip();
        self.u_vel.flip();
        self.v_vel.flip();
    }

    pub fn add_inflow(&self, x: f64, y: f64, width: f64, height: f64, d: f64, u: f64, v: f64) {
        self.density.add_inflow();
        self.u_vel.add_inflow();
        self.v_vel.add_inflow();
    }

    pub fn max_timestep(&self) -> f64 {
        let mut max_velocity: f64 = 0.0;

        for y in 0..self.height {
            for x in 0..self.width {
                let u: f64 = self.u_vel.lerp(x as f64 + 0.5, y as f64 + 0.5);
                let v: f64 = self.v_vel.lerp(x as f64 + 0.5, y as f64 + 0.5);

                let velocity: f64 = (u * u + v * v).sqrt();
                if max_velocity < velocity {
                    max_velocity = velocity
                }
            }
        }

        let mut max_timestep = 2.0 * self.hx / max_velocity;

        if max_timestep > 1.0 {
            max_timestep = 1.0;
        }
        max_timestep
    }

    pub fn to_image(&self, buffer: &mut Box<[u8; SIZE_X * SIZE_Y * 4]>) {
        for i in 0..(self.width * self.height) {
            let mut shade: u8 = ((1.0 - self.density.src[i]) * 255.0) as u8;
            shade = max(min(shade, 255), 0);

            buffer[i * 4 + 0] = shade;
            buffer[i * 4 + 1] = shade;
            buffer[i * 4 + 2] = shade;
            buffer[i * 4 + 3] = 0xFF;
        }
    }
}


const SIZE_X: usize = 128;
const SIZE_Y: usize = 128;

fn main() {
    const DENSITY: f64 = 0.1;
    const TIMESTEP: f64 = 0.005;

    let mut image = Box::new([0u8; SIZE_X * SIZE_Y * 4]);

    let mut solver = FluidSolver::new(SIZE_X, SIZE_Y, DENSITY);

    let mut time: f64 = 0.0;
    let mut iterations = 0;

    while time < 1.0 {
        solver.add_inflow(0.45, 0.2, 0.1, 0.01, 1.0, 0.0, 3.0);
        solver.update(TIMESTEP);
        time += TIMESTEP;
        stdout().flush();
        solver.to_image(&mut image);
        println!("Encoding {}.", iterations);
        iterations += 1;

        image::save_buffer(format!("output/{}.png", iterations), &image.to_vec(), SIZE_X as u32, SIZE_Y as u32, image::RGBA(8));
    }



}