#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>

const double MU = 3.986004418e14; // m^3/s^2
const double R_EARTH = 6378137.0;  // m
const double PI = 3.14159265358979323846;

const double J2 = 1.08263e-3;
const double CD = 2.2;
const double AREA = 1.0;
const double RHO_400KM = 4.0e-12;
const double SCALE_HEIGHT = 50000.0;

struct Thruster {
    double thrust_N = 0.050; // 50 mN
    double isp = 1500.0;     // s
    double g0 = 9.80665;

    double m_dot() const {
        return thrust_N / (isp * g0);
    }
};

// 13D State Vector: [x, y, z, vx, vy, vz, m, lx, ly, lz, lvx, lvy, lvz]
// Includes 6 Costate variables (lx...lvz) for optimal thrust vectoring
struct FullState {
    double x, y, z;
    double vx, vy, vz;
    double m;
    double lx, ly, lz;       // Position costates
    double lvx, lvy, lvz;   // Velocity costates (Lawden's Primer Vector)

    FullState operator+(const FullState& b) const {
        return {x + b.x, y + b.y, z + b.z, 
                vx + b.vx, vy + b.vy, vz + b.vz, 
                m + b.m,
                lx + b.lx, ly + b.ly, lz + b.lz,
                lvx + b.lvx, lvy + b.lvy, lvz + b.lvz};
    }

    FullState operator*(double scalar) const {
        return {x * scalar, y * scalar, z * scalar, 
                vx * scalar, vy * scalar, vz * scalar, 
                m * scalar,
                lx * scalar, ly * scalar, lz * scalar,
                lvx * scalar, lvy * scalar, lvz * scalar};
    }
};

void cross_product(const double a[3], const double b[3], double out[3]) {
    out[0] = a[1] * b[2] - a[2] * b[1];
    out[1] = a[2] * b[0] - a[0] * b[2];
    out[2] = a[0] * b[1] - a[1] * b[0];
}

double calculate_inclination_deg(const FullState& s) {
    double r[3] = {s.x, s.y, s.z};
    double v[3] = {s.vx, s.vy, s.vz};
    double h[3];
    cross_product(r, v, h);

    double h_mag = std::sqrt(h[0]*h[0] + h[1]*h[1] + h[2]*h[2]);
    return std::acos(h[2] / h_mag) * (180.0 / PI);
}

double calculate_density(double altitude_m) {
    if (altitude_m > 1000000.0) return 0.0;
    return RHO_400KM * std::exp(-(altitude_m - 400000.0) / SCALE_HEIGHT);
}

// Compute derivatives for physical motion AND costate dynamics
FullState get_derivatives_optimal(const FullState& s, const Thruster& engine) {
    double r_mag = std::sqrt(s.x * s.x + s.y * s.y + s.z * s.z);
    double v_mag = std::sqrt(s.vx * s.vx + s.vy * s.vy + s.vz * s.vz);
    double altitude = r_mag - R_EARTH;

    // 1. Gravity + J2 + Drag Accelerations
    double ax_g = -MU * s.x / (r_mag * r_mag * r_mag);
    double ay_g = -MU * s.y / (r_mag * r_mag * r_mag);
    double az_g = -MU * s.z / (r_mag * r_mag * r_mag);

    double z_sq_over_r_sq = (s.z * s.z) / (r_mag * r_mag);
    double j2_factor = (1.5 * J2 * MU * R_EARTH * R_EARTH) / std::pow(r_mag, 5);

    double ax_j2 = j2_factor * s.x * (5.0 * z_sq_over_r_sq - 1.0);
    double ay_j2 = j2_factor * s.y * (5.0 * z_sq_over_r_sq - 1.0);
    double az_j2 = j2_factor * s.z * (5.0 * z_sq_over_r_sq - 3.0);

    double rho = calculate_density(altitude);
    double drag_factor = -0.5 * rho * (CD * AREA / s.m) * v_mag;
    double ax_drag = drag_factor * s.vx;
    double ay_drag = drag_factor * s.vy;
    double az_drag = drag_factor * s.vz;

    // 2. Pontryagin's Optimal Steering Vector from Velocity Costates (Primer Vector)
    double lv_mag = std::sqrt(s.lvx * s.lvx + s.lvy * s.lvy + s.lvz * s.lvz);
    
    // Optimal thrust direction is anti-parallel to velocity costate
    double ux = -s.lvx / lv_mag;
    double uy = -s.lvy / lv_mag;
    double uz = -s.lvz / lv_mag;

    double accel_thrust = engine.thrust_N / s.m;
    double ax_t = accel_thrust * ux;
    double ay_t = accel_thrust * uy;
    double az_t = accel_thrust * uz;

    // 3. Physical State Derivatives
    double dx = s.vx;
    double dy = s.vy;
    double dz = s.vz;
    double dvx = ax_g + ax_j2 + ax_drag + ax_t;
    double dvy = ay_g + ay_j2 + ay_drag + ay_t;
    double dvz = az_g + az_j2 + az_drag + az_t;
    double dm = -engine.m_dot();

    // 4. Costate Differential Equations (Adjoint Dynamics)
    // dl_r / dt = - dH/dr  |  dl_v / dt = - dH/dv
    double dlx = (MU / std::pow(r_mag, 3)) * s.lvx - (3.0 * MU * s.x * (s.x * s.lvx + s.y * s.lvy + s.z * s.lvz)) / std::pow(r_mag, 5);
    double dly = (MU / std::pow(r_mag, 3)) * s.lvy - (3.0 * MU * s.y * (s.x * s.lvx + s.y * s.lvy + s.z * s.lvz)) / std::pow(r_mag, 5);
    double dlz = (MU / std::pow(r_mag, 3)) * s.lvz - (3.0 * MU * s.z * (s.x * s.lvx + s.y * s.lvy + s.z * s.lvz)) / std::pow(r_mag, 5);

    double dlvx = -s.lx;
    double dlvy = -s.ly;
    double dlvz = -s.lz;

    return {dx, dy, dz, dvx, dvy, dvz, dm, dlx, dly, dlz, dlvx, dlvy, dlvz};
}

FullState rk4_step_optimal(const FullState& s, const Thruster& engine, double dt) {
    FullState k1 = get_derivatives_optimal(s, engine);
    FullState k2 = get_derivatives_optimal(s + k1 * (dt / 2.0), engine);
    FullState k3 = get_derivatives_optimal(s + k2 * (dt / 2.0), engine);
    FullState k4 = get_derivatives_optimal(s + k3 * dt, engine);

    return s + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (dt / 6.0);
}

int main() {
    double initial_alt_km = 400.0;
    double target_alt_km = 1000.0;
    
    double r0 = R_EARTH + (initial_alt_km * 1000.0);
    double v0 = std::sqrt(MU / r0);

    // Initial state with optimal costate initialization for prograde spiral
    FullState state = {
        r0, 0.0, 0.0,           // Position (x, y, z)
        0.0, v0, 0.0,           // Velocity (vx, vy, vz)
        100.0,                  // Mass (kg)
        0.0, 0.0, 0.0,          // Position Costates (lx, ly, lz)
        0.0, -1.0, -0.1         // Velocity Costates (lvx, lvy, lvz) -> steering bias
    }; 

    Thruster engine;

    double dt = 10.0;
    double current_time = 0.0;
    double target_radius = R_EARTH + (target_alt_km * 1000.0);

    std::ofstream data_file("optimal_trajectory.csv");
    data_file << "time_days,x_km,y_km,z_km,altitude_km,inclination_deg,mass_kg\n";

    std::cout << "Starting Optimal Control (Costate / Primer Vector) Simulation..." << std::endl;

    while (true) {
        double current_r = std::sqrt(state.x * state.x + state.y * state.y + state.z * state.z);
        
        if (current_r >= target_radius || state.m <= 0.0) {
            break;
        }

        state = rk4_step_optimal(state, engine, dt);
        current_time += dt;

        if (static_cast<int>(current_time) % 3600 == 0) {
            double time_days = current_time / 86400.0;
            double altitude_km = (current_r - R_EARTH) / 1000.0;
            double inc_deg = calculate_inclination_deg(state);

            data_file << time_days << "," << state.x/1000.0 << "," 
                      << state.y/1000.0 << "," << state.z/1000.0 << "," 
                      << altitude_km << "," << inc_deg << "," << state.m << "\n";
        }
    }

    data_file.close();

    std::cout << "--- OPTIMAL TRAJECTORY SIMULATION COMPLETE ---" << std::endl;
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Transfer Time:       " << (current_time / 86400.0) << " days" << std::endl;
    std::cout << "Final Inclination:   " << calculate_inclination_deg(state) << " deg" << std::endl;
    std::cout << "Propellant Used:     " << (100.0 - state.m) << " kg" << std::endl;
    std::cout << "Data saved to 'optimal_trajectory.csv'" << std::endl;

    return 0;
}