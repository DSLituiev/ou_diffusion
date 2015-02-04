#include <iostream>
#include <string.h> // memset
#include <stdio.h>
#include <iomanip>      // std::setprecision
#include <sqlite3pp.h>
#include <math.h>       /* sqrt */

using namespace std;

class Params{
public:
    int const N ;
    int const T ;
    const size_t GRID_X_POINTS ;

    size_t dimensionality = 1;
    float D ;
    float dt;
    float k;
    int const RANDOMNESS_BITS_U = 10;
    float kT = 1.3806488131313e-23 * 293.15 * 1e24;
    float friction;

    float VARIANCE, sigma_g;
    float sqrt_dt;
    float const GRID_DX;
 
//    const char * RUN_NAME;
//    const char * OUT_FILE ;
    string OUT_FILE_STR ;
    string RUN_NAME_STR ;

    Params(int n, int t, int gp, float d, float k_, float dt_, float h_, const char * rn, const char * out): \
    N(n), T(t), GRID_X_POINTS(gp), D(d), k(k_), dt(dt_), GRID_DX(h_), RUN_NAME_STR(rn), OUT_FILE_STR(out)
    {
        // Without friction
        VARIANCE = 2.0 * (float) dimensionality * D;
        sigma_g = sqrt(VARIANCE);
/*
        sigma_W = sqrt(dt);
        sigma_U = 1/3*dt*dt*dt;
        rho_W_U = 1/2*dt*dt;
*/
        friction = D / kT;
        sqrt_dt = sqrt(dt);

        clog << "database : " << OUT_FILE_STR << endl;
        printf("D  : %1.2g\n" , D) ; //clog << "D  : " << D << endl;
        printf("friction  : %1.2g\n" , friction) ; //clog << "D  : " << D << endl;
        clog << "k  : " << k << endl;
        clog << "dt : " << dt << endl; 
        clog << "variance : " << VARIANCE << endl; 
    }

    Params(int n, int t, int gp, float d, float k_, float dt_, float h_, float f_, const char * rn, const char * out): \
    N(n), T(t), GRID_X_POINTS(gp), D(d), k(k_), dt(dt_), GRID_DX(h_), RUN_NAME_STR(rn), OUT_FILE_STR(out), friction(f_)
    {
        // Without friction
        VARIANCE = 2.0 * (float) dimensionality * D;
        sigma_g = sqrt(VARIANCE);
/*
        sigma_W = sqrt(dt);
        sigma_U = 1/3*dt*dt*dt;
        rho_W_U = 1/2*dt*dt;
*/
        sqrt_dt = sqrt(dt);

        clog << "database : " << OUT_FILE_STR << endl;
        printf("D  : %1.2g\n" , D) ; //clog << "D  : " << D << endl;
        printf("friction  : %1.2g\n" , friction) ; //clog << "D  : " << D << endl;
        clog << "k  : " << k << endl;
        clog << "dt : " << dt << endl; 
        clog << "variance : " << VARIANCE << endl; 
    }

    ~Params(){
    }


void register_table(sqlite3pp::database * db) const {
    bool drop = false;
    if (drop){
        std::string sql_drop( "DROP TABLE IF EXISTS register");
        db->execute(sql_drop.c_str());
    }

    char sql[] = "CREATE TABLE IF NOT EXISTS register ("  \
              "id  TEXT PRIMARY KEY, " \
              "particles INT, " \
              "grid_points INT, " \
              "time_points INT, " \
              "diffusion FLOAT , " \
              "young FLOAT, " \
              "dt FLOAT, "  \
              "grid_dx FLOAT " \
              ");";
    try {
         if ( db->execute(sql) != 0) {
            clog << "[sqlite3 failed:] " << sql << endl ;
            cerr << db->error_msg() << endl;
         }
    } catch (exception& ex) {
         cerr << ex.what() << endl;
    }
    
    char sql_ins[] = "INSERT OR REPLACE INTO register "  \
                  "(id, particles, grid_points, time_points, diffusion, young, dt, grid_dx) " \
           "VALUES (:id, :particles, :grid_points, :time_points, :diffusion, :young, :dt, :dx); ";

    sqlite3pp::command cmd( *db, sql_ins);
    cmd.bind(":id", RUN_NAME_STR.c_str() );
    cmd.bind(":particles", N );
    cmd.bind(":grid_points", (int) GRID_X_POINTS );
    cmd.bind(":time_points", T );
    cmd.bind(":diffusion", D );
    cmd.bind(":young", k );
    cmd.bind(":dt", dt );
    cmd.bind(":dx", GRID_DX );
    try {
         if ( cmd.execute() != 0) {
            clog << "[sqlite3 failed:] " << sql_ins << endl ;
             clog << RUN_NAME_STR.c_str() << "\t" ; 
             cerr << db->error_msg() << endl;
         }
    } catch (exception& ex) {
         cerr << ex.what() << endl;
    }
  }

    float drift_term(const float & x, const float & y) const{
        return  y - k * x * dt ;
    }

    float drift_term_velocity(const float & x, const float &y) const{
        return  -k * x  - friction* y;
    }

    float diffusion_term(const float & rndn) const{
        return rndn * sqrt_dt * sigma_g;
    }


    int map_x_to_grid(float & x) const {
        int x_ = (int) round(x / GRID_DX);
        int gg;
        gg = (abs(x_) < GRID_X_POINTS / 2) ? ( GRID_X_POINTS / 2 - x_) : (x_<0 ? 0: GRID_X_POINTS-1);
        return gg;
    }    

};
