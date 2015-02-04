#include <iostream>
#include <string.h> // memset
#include <stdio.h>
#include <stdlib.h>     /* abs */
#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <iomanip>      // std::setprecision
#include <sqlite3pp.h>
#include "tclap/CmdLine.h"
#include "Params.h"

// #include "gsl_randist.h"

using namespace std;
using namespace TCLAP;

void uint16_to_char(char * out_buffer, const unsigned short int & val){
    out_buffer[0] = (val) & 0xFF;  // low byte
    out_buffer[1] = (val >> 8) & 0xFF;  // high byte
    return;
}

void array_uint16_to_char(char * out_buffer, unsigned short int * arr, size_t len){
    for (size_t ii = 0; ii<len; ii++){
        uint16_to_char( out_buffer + 2*ii, arr[ii]);
    }
    return;
}

// http://stackoverflow.com/questions/7217791/random-numbers-in-c0x
/*
int map_x_to_grid(float & x, size_t const & GRID_X_POINTS){
    int x_ = (int) round(x);
    int gg;
    gg = (abs(x_) < GRID_X_POINTS / 2) ? ( GRID_X_POINTS / 2 - x_) : (x_<0 ? 0: GRID_X_POINTS-1);
    return gg;
}

int map_x_to_grid(float & x, size_t const & GRID_X_POINTS, float const & GRID_DX){
    int x_ = (int) round(x / GRID_DX);
    int gg;
    gg = (abs(x_) < GRID_X_POINTS / 2) ? ( GRID_X_POINTS / 2 - x_) : (x_<0 ? 0: GRID_X_POINTS-1);
    return gg;
}
*/
/*
void elastic_acceleration( float dxdt[], float k, float * positions, size_t N ){
   
    for (size_t nn=0; nn <N; nn++){
        dxdt[nn] = k * positions[nn];
    }
    return;
} */
 
void init_sql_table(sqlite3pp::database * db, const char run_name [], bool drop){
    if (drop){
        std::string sql_drop( "DROP TABLE IF EXISTS ");
        std::string name( run_name );
        sql_drop = sql_drop + run_name;
        db->execute(sql_drop.c_str());
    }

    char sql[256];
    sprintf(sql, "CREATE TABLE IF NOT EXISTS %s ("  \
              "time INT PRIMARY KEY     NOT NULL , " \
              "grid BLOB, " \
              "track1 FLOAT, " \
              "track2 FLOAT, " \
              "track3 FLOAT );", run_name );
    try {
         if (db->execute(sql)!=0){
            clog << "[sqlite3 failed:] " << sql << endl;
            cerr << db->error_msg() << endl;
         }
    } catch (exception& ex) {
         cerr << ex.what() << endl;
    }
        // clog << "table " << id+1 << " has been created" << endl;

}

void print_to_db(sqlite3pp::database * db, const char * run_name, \
     unsigned short int * grid_t, float *tracks, size_t & tt, const size_t & GRID_X_POINTS ){

    char buffer[GRID_X_POINTS*2];
    bool PRINT_FLAG = false;    
    array_uint16_to_char(buffer, grid_t, GRID_X_POINTS);

    char sql[256];
    sprintf (sql, "INSERT INTO %s "  \
                  "(time, grid, track1, track2, track3) " \
                  "VALUES (:time, :grid, :t1, :t2, :t3) ", \
                   run_name );

    sqlite3pp::command cmd( *db, sql);
    cmd.bind(":time", (int) tt);
//  cmd.bind(":tot_cov", buffer );
//    int blob_status = cmd.bind(":tot_cov", (void const*) buffer, (int) GRID_X_POINTS*2 );
    int blob_status = cmd.bind(":grid", (void const*) grid_t, (int) GRID_X_POINTS*2 );
    cmd.bind(":t1", tracks[0] );
    cmd.bind(":t2", tracks[1] );
    cmd.bind(":t3", tracks[2] );
    try {
        cmd.execute();
    } catch (exception& ex) {
        cerr << ex.what() << endl;
    }

//    std::cout << "printed time point " << tt << " points: " <<  GRID_X_POINTS*2 <<" status :" << blob_status << std::endl;
  if (PRINT_FLAG ){
    std::cout << "time point: " << tt << std::endl;
    for (size_t gg=0; gg<GRID_X_POINTS ; gg++){
        std::cout << setfill(' ') << setw(3) << grid_t[gg] << "  ";
    }
    std::cout << std::endl;
  }
}

Params * parse_cmd_line(int argc, char *argv[]) {

    CmdLine cmdss("Comp Motiv", ' ', "0");

    ValueArg<string> out("o", "out",
            "database path", false, "#", "string");
    cmdss.add(out);

    ValueArg<int> N_("n", "particles", "number of particles", true,
            0, "int");
    cmdss.add(N_);

    ValueArg<string> l_("l", "label", "run / table label", false,
            "simple", "float");
    cmdss.add(l_);

    ValueArg<int> T_("t", "time", "time points", true,
            0, "int");
    cmdss.add(T_);

    ValueArg<int> G_("g", "grid", "grid points", true,
            0, "int");
    cmdss.add(G_);

    ValueArg<float> D_("d", "diffusion", "diffusion coefficient", true,
            0.0, "float");
    cmdss.add(D_);

    ValueArg<float> dt_("s", "step", "step size", true,
            1.0, "float");
    cmdss.add(dt_);

    ValueArg<float> k_("k", "young", "Young modulus", false,
            0.0, "float");
    cmdss.add(k_);

    ValueArg<float> f_("f", "friction", "friction coefficient", false,
            0.0, "float");
    cmdss.add(f_);
    
    ValueArg<float> r_("r", "grid_dx", "grid spacing", false,
            1.0, "float");
    cmdss.add(r_);


    Params * p;
    
    try {
        cmdss.parse(argc, argv); //parse arguments
        std::string output = out.getValue();
        if (output[0] == '#') {
            output = "plots.db";
        }
        int N = N_.getValue();
        int T = T_.getValue();
        int G = G_.getValue();

        if (f_.getValue()>0){
             p = new Params(N, T, G, D_.getValue(), k_.getValue(), dt_.getValue(), r_.getValue(), f_.getValue(), \
             l_.getValue().c_str(), output.c_str() );
        }else{
            p = new Params(N, T, G, D_.getValue(), k_.getValue(), dt_.getValue(), r_.getValue(), l_.getValue().c_str(), output.c_str() );
        }
    } catch (ArgException &e) {
        StdOutput out;
        out.failure(cmdss, e);
    }
    return p;
}

int main(int argc, char *argv[]) {
    // char OUT_FILE [] = "plots.db";
    // char RUN_NAME[] = "elastic";
    std::cout << "starting... ";
   
    const Params * p = parse_cmd_line(argc, &*argv);
   
    // cout << "variance : " << p->VARIANCE << endl;
   
    // Params p = Params(2000, 5000, 400, 5.0, 0.0, 0.1);
    // params p = params(N, T, GRID_X_POINTS, D, k, dt);

    // DATABASE INITIALISATION
    sqlite3pp::database * db;
    cout << "database : " << p->OUT_FILE_STR.c_str() << endl;
    cout << "table : " << p->RUN_NAME_STR.c_str() << endl;
    db  = new sqlite3pp::database( p->OUT_FILE_STR.c_str() ) ;

    sqlite3pp::transaction * xct = new sqlite3pp::transaction(*db);

    p->register_table(db);
    init_sql_table(db, p->RUN_NAME_STR.c_str(), true);

//
    float ** particles = new float* [p->T];
    for (size_t tt=0; tt< p->T; tt++){
        particles[tt] = new float [p->N];
        memset(particles[tt], 0, sizeof(float) * (p->N) );
    }

    float * velocity_particles = new float [p->N];
    memset(velocity_particles, 0, sizeof(float) * (p->N) );

    float * dv = new float [p->N];

    unsigned short int ** grid = new unsigned short int*[p->T];
    
    for (size_t tt=0; tt<p->T; tt++){
        grid[tt] = new unsigned short int [p->GRID_X_POINTS];
        memset(grid[tt], 0, sizeof(unsigned short int) * p->GRID_X_POINTS);
    }

    for (size_t nn=0; nn<p->N ; nn++){
        grid[0 ][ p->map_x_to_grid(particles[0][nn]) ]++;
    }
    
/*    for (size_t gg=0; gg < p->GRID_X_POINTS ; gg++){
        std::cout << grid[0][gg] << "  ";
    }
    cout << endl;
*/
  
    size_t tt = 0;
    print_to_db( db, p->RUN_NAME_STR.c_str(), grid[tt], particles[tt], tt, p->GRID_X_POINTS);
/*
    for (size_t gg=0; gg<GRID_X_POINTS ; gg++){
        std::cout << setfill(' ') << setw(3) << grid[0][gg] << "  ";
    }
    cout << endl;
    cout << "===========================================";
    cout << endl;
*/
    // Good random seed, good engine
    auto rnd1 = std::mt19937(std::random_device{}());
//    auto rnd2 = std::mt19937(std::random_device{}());
/*    const gsl_rng_type * gslT;
    gsl_rng * gslr;
    gsl_rng_env_setup();
    gslT = gsl_rng_default;
    gslr = gsl_rng_alloc(gslT);
*/    
    float theta;
    float dW;
    float F1y, F2y; // Honeycutt 1992

    std::normal_distribution<> normal_dist(0, 1);
    
    std::cout << "initialization completed." << std::endl;
    
    cout << grid[tt][0] << "\t" ;

    float dvdx = - p->k * p->dt;

    for (size_t tt=1; tt < p->T; tt++){
        for (size_t nn=0; nn < p->N ; nn++){
            // theta = std::generate_canonical<double, p->RANDOMNESS_BITS_U>(rnd2) ;

            dW = normal_dist(rnd1);
             
        //    gsl_ran_bivariate_gaussian (gslr, double sigma_W, double sigma_U, double rho, &dW, &dU)

           // elastic_acceleration(dv, k, particles[tt-1], N);

//            dv[nn] = p->drift_term(particles[tt-1][nn]);

            
            // EXPLICIT EULER
            //particles[tt][nn] = particles[tt-1][nn] + p->diffusion_term(dW) + velocity_particles[nn]  * p->dt;
            //velocity_particles[nn] +=  p->drift_term_velocity(particles[tt-1][nn]) * p->dt;
 
            // Milstein 
            //  p->sqrt_dt * p->dt * dW * dvdx * 1/2 + \
            //  p->dt * p->dt * 1/2* dvdx *  velocity_particles[nn] ;

            // Euler - Heun
            // particles[tt][nn] = particles[tt-1][nn] + dW * p->sqrt_dt  + velocity_particles[nn] * p->dt;

            // Runge-Kutta-Honeycutt
            
            particles[tt][nn] = particles[tt-1][nn] + p->dt * velocity_particles[nn] ; 

            velocity_particles[nn] += p->dt * p->drift_term_velocity(particles[tt-1][nn], velocity_particles[nn]) \
                                    + p->diffusion_term( dW );
            
//             int grid_point = p->map_x_to_grid( particles[tt][nn] );

             grid[tt][ p->map_x_to_grid( particles[tt][nn] ) ]++; 
        };
//        cout << ' ' <<  setfill(' ') << setw(5)  << tt;
        printf( " %1.3f | ", particles[tt][20]);
        if (tt % 10 == 0){ cout << endl;}

        print_to_db( db, p->RUN_NAME_STR.c_str(), grid[tt], particles[tt], tt , p->GRID_X_POINTS);
//        if (tt % 100 == 0){
//            cout << '\r' <<  setfill(' ') << setw(5)  << tt;
//        }
    }

    //gsl_rng_free (gslr); // GSL random number generator

    xct->commit(); 
    
    std::cout << std::endl;
    std::cout << "finished successfully";
    std::cout << std::endl;

    for (size_t tt=0; tt<p->T; tt++){
        delete [] particles[tt];
        delete [] grid[tt];
    }
    delete [] grid;
    delete [] particles;


    return 0;
}


