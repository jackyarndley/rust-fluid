use crate::util::rotate;

enum CellType {
    FLUID,
    SOLID
}

struct SolidBody {
    pos_x: f64,
    pos_y: f64,
    scale_x: f64,
    scale_y: f64,
    theta: f64,
    vel_x: f64,
    vel_y: f64,
    omega: f64
}

impl SolidBody {
    pub fn new(pos_x: f64, pos_y: f64, scale_x: f64, scale_y: f64, theta: f64, vel_x: f64, vel_y: f64, omega: f64) -> SolidBody {
        SolidBody {
            pos_x,
            pos_y,
            scale_x,
            scale_y,
            theta,
            vel_x,
            vel_y,
            omega
        }
    }

    fn global_to_local(&self, x: &mut f64, y: &mut f64) {
        x -= self.pos_x;
        y -= self.pos_y;
        rotate(x, y, -self.theta);
        x /= self.scale_x;
        y /= self.scale_y;
    }

    fn local_to_global(&self, x: &mut f64, y: &mut f64) {
        x *= self.scale_x;
        y *= self.scale_y;
        rotate(x, y, self.theta);
        x += self.pos_x;
        y += self.pos_y;
    }
}

