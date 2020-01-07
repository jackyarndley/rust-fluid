use crate::util::helper::{length, max, rotate, nsgn};

pub enum SolidType {
    Box,
    Sphere
}

pub struct SolidBody {
    pos_x: f64,
    pos_y: f64,
    scale_x: f64,
    scale_y: f64,
    theta: f64,
    vel_x: f64,
    vel_y: f64,
    vel_theta: f64,
    solid_type: SolidType
}

impl SolidBody {
    pub fn new_box(pos_x: f64, pos_y: f64, scale_x: f64, scale_y: f64, theta: f64, vel_x: f64, vel_y: f64, vel_theta: f64) -> Self {
        SolidBody {
            pos_x,
            pos_y,
            scale_x,
            scale_y,
            theta,
            vel_x,
            vel_y,
            vel_theta,
            solid_type: SolidType::Box
        }
    }

    pub fn new_sphere(pos_x: f64, pos_y: f64, scale: f64, theta: f64, vel_x: f64, vel_y: f64, vel_theta: f64) -> Self {
        SolidBody {
            pos_x,
            pos_y,
            scale_x: scale,
            scale_y: scale,
            theta,
            vel_x,
            vel_y,
            vel_theta,
            solid_type: SolidType::Sphere
        }
    }

    fn global_to_local(&self, x: &mut f64, y: &mut f64) {
        *x -= self.pos_x;
        *y -= self.pos_y;
        rotate(x, y, -self.theta);
        *x /= self.scale_x;
        *y /= self.scale_y;
    }

    fn local_to_global(&self, x: &mut f64, y: &mut f64) {
        *x *= self.scale_x;
        *y *= self.scale_y;
        rotate(x, y, self.theta);
        *x += self.pos_x;
        *y += self.pos_y;
    }

    pub fn velocity_x(&self, _x: f64, y: f64) -> f64 {
        (self.pos_y - y) * self.vel_theta + self.vel_x
    }

    pub fn velocity_y(&self, x: f64, _y: f64) -> f64 {
        (self.pos_x - x) * self.vel_theta + self.vel_y
    }

    pub fn velocity(&self, vx: &mut f64, vy: &mut f64, x: f64, y: f64) {
        *vx = self.velocity_x(x, y);
        *vy = self.velocity_y(x, y);
    }

    pub fn update(&mut self, timestep: f64) {
        self.pos_x += self.vel_x * timestep;
        self.pos_y += self.vel_y * timestep;
        self.theta += self.vel_theta * timestep;
    }

    pub fn distance(&self, mut x: f64, mut y: f64) -> f64 {
        match self.solid_type {
            SolidType::Box => {
                x -= self.pos_x;
                y -= self.pos_y;
                rotate(&mut x, &mut y, -self.theta);
                let dx = x.abs() - self.scale_x * 0.5;
                let dy = y.abs() - self.scale_y * 0.5;

                if (dx >= 0.0) || (dy >= 0.0) {
                    length(max(dx, 0.0), max(dy, 0.0))
                } else {
                    max(dx, dy)
                }
            },
            SolidType::Sphere => {
                length(x - self.pos_x, y - self.pos_y)
            }
        }
    }

    pub fn closest_surface_point(&self, x: &mut f64, y: &mut f64) {
        match self.solid_type {
            SolidType::Box => {
                *x -= self.pos_x;
                *y -= self.pos_y;
                rotate(x, y, -self.theta);
                let dx = x.abs() - self.scale_x * 0.5;
                let dy = y.abs() - self.scale_y * 0.5;

                if dx > dy {
                    *x = nsgn(*x) * 0.5 * self.scale_x;
                } else {
                    *y = nsgn(*y) * 0.5 * self.scale_y;
                }

                rotate(x, y, self.theta);
                *x += self.pos_x;
                *y += self.pos_y;
            },
            SolidType::Sphere => {
                self.global_to_local(x, y);
                let r = length(*x, *y);

                if r < 1e-4 {
                    *x = 0.5;
                    *y = 0.0;
                } else {
                    *x /= 2.0 * r;
                    *y /= 2.0 * r
                }

                self.local_to_global(x, y);
            }
        }
    }

    pub fn distance_normal(&self, mut x: f64, mut y: f64) -> (f64, f64) {
        let mut normal_x: f64;
        let mut normal_y: f64;
        match self.solid_type {
            SolidType::Box => {
                x -= self.pos_x;
                y -= self.pos_y;
                rotate(&mut x, &mut y, -self.theta);

                if x.abs() - self.scale_x * 0.5 > y.abs() - self.scale_y * 0.5 {
                    normal_x = nsgn(x);
                    normal_y = 0.0;
                } else {
                    normal_x = 0.0;
                    normal_y = nsgn(y);
                }

                rotate(&mut normal_x, &mut normal_y, self.theta);

                (normal_x, normal_y)
            },
            SolidType::Sphere => {
                x -= self.pos_x;
                y -= self.pos_y;
                let r = length(x, y);

                if r < 1e-4 {
                    normal_x = 1.0;
                    normal_y = 0.0;
                } else {
                    normal_x = x/r;
                    normal_y = y/r;
                }

                (normal_x, normal_y)
            }
        }
    }
}

