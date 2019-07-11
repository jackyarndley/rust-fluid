use crate::util::{rotate, max, length};

enum CellType {
    FLUID,
    SOLID
}

enum SolidType {
    BOX,
    SPHERE
}

struct SolidBody {
    pos_x: f64,
    pos_y: f64,
    scale_x: f64,
    scale_y: f64,
    theta: f64,
    vel_x: f64,
    vel_y: f64,
    omega: f64,
    solid_type: SolidType
}

impl SolidBody {
    pub fn new(pos_x: f64, pos_y: f64, scale_x: f64, scale_y: f64, theta: f64, vel_x: f64, vel_y: f64, omega: f64, solid_type: SolidType) -> Self {
        Self {
            pos_x,
            pos_y,
            scale_x,
            scale_y,
            theta,
            vel_x,
            vel_y,
            omega,
            solid_type
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

    fn distance(&self, mut x: f64, mut y: f64) -> f64 {
        match self.solid_type {
            SolidType::BOX => {
                x -= self.pos_x;
                y -= self.pos_y;

                rotate(&mut x, &mut y, -self.theta);

                let dx = x.abs() - 0.5 * self.scale_x;
                let dy = y.abs() - 0.5 * self.scale_y;

                if dx >= 0.0 || dy >= 0.0 {
                    length(max(dx, 0.0), max(dy, 0.0))
                } else {
                    max(dx, dy)
                }
            },
            SolidType::SPHERE => {
                length(x - self.pos_x, y - self.pos_y) - self.scale_x * 0.5
            }
        }
    }

    fn closest_surface_point(&self, x: &mut f64, y: &mut f64) {
        match self.solid_type {
            SolidType::BOX => {
                *x -= self.pos_x;
                *y -= self.pos_y;

                rotate(x, y, -self.theta);

                let dx = x.abs() - 0.5 * self.scale_x;
                let dy = y.abs() - 0.5 * self.scale_y;

                if dx > dy {
                    *x = 0.5 * x.signum() * self.scale_x;
                } else {
                    *y = 0.5 * y.signum() * self.scale_y;
                }

                rotate(x, y, self.theta);

                *x += self.pos_x;
                *y += self.pos_y;
            },
            SolidType::SPHERE => {
                self.global_to_local(x, y);

                let mut r: f64 = length(*x, *y);

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

    fn distance_normal(&self, nx: &mut f64, ny: &mut f64, mut x: f64, mut y: f64) {
        match self.solid_type {
            SolidType::BOX => {
                x -= self.pos_x;
                y -= self.pos_y;

                rotate(&mut x, &mut y, -self.theta);

                if x.abs() - 0.5 * self.scale_x > y.abs() - 0.5 * self.scale_y {
                    *nx = x.signum();
                    *ny = 0.0;
                } else {
                    *nx = 0.0;
                    *ny = y.signum();
                }

                rotate(nx, ny, self.theta)
            },
            SolidType::SPHERE => {
                self.global_to_local(&mut x, &mut y);

                let mut r: f64 = length(x, y);

                if r < 1e-4 {
                    *nx = 1.0;
                    *ny = 0.0;
                } else {
                    *nx = x / r;
                    *ny = y / r
                }

                self.local_to_global(&mut x, &mut y)
            }
        }
    }

}