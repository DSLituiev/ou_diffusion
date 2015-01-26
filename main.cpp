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

using namespace std;
using namespace TCLAP;
// size_t GRID_X_POINTS;


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
int map_x_to_grid(float & x, size_t const & GRID_X_POINTS){
    int x_ = (int) round(x);
    int gg;
    gg = (abs(x_) < GRID_X_POINTS / 2) ? ( GRID_X_POINTS / 2 - x_) : (x_<0 ? 0: GRID_X_POINTS-1);
    return gg;
    }

void elastic_acceleration( float dxdt[], float k, float * positions, size_t N ){
   
    for (size_t nn=0; nn <N; nn++){
        dxdt[nn] = k * positions[nn];
    }
    return;
}

 
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
              "grid BLOB);", run_name );
    try {
         clog << "[sqlite3 :] " << sql << " >>> " ;
         clog << db->execute(sql) << endl;
    } catch (exception& ex) {
         cerr << ex.what() << endl;
    }
        // clog << "table " << id+1 << " has been created" << endl;

}

void print_to_db(sqlite3pp::database * db, const char * run_name, unsigned short int * grid_t, size_t & tt, const size_t & GRID_X_POINTS ){
    char buffer[GRID_X_POINTS*2];
    bool PRINT_FLAG = false;    
    array_uint16_to_char(buffer, grid_t, GRID_X_POINTS);

    char sql[256];
    sprintf (sql, "INSERT INTO %s "  \
                  "(time, grid) " \
                  "VALUES (:time, :grid) ", \
                   run_name );

    sqlite3pp::command cmd( *db, sql);
    cmd.bind(":time", (int) tt);
//  cmd.bind(":tot_cov", buffer );
//    int blob_status = cmd.bind(":tot_cov", (void const*) buffer, (int) GRID_X_POINTS*2 );
    int blob_status = cmd.bind(":grid", (void const*) grid_t, (int) GRID_X_POINTS*2 );

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
        int D = D_.getValue();

        p = new Params(N, T, G, D, k_.getValue(), dt_.getValue(), l_.getValue().c_str(), output.c_str() );

    } catch (ArgException &e) {
        StdOutput out;
        out.failure(cmdss, e);
    }
return p;
}
/*
void diffusor(){
    for (size_t nn=0; nn < p->N ; nn++){
        // theta = std::generate_canonical<double, p->RANDOMNESS_BITS_U>(rnd2) ;

        dr = normal_dist(rnd1);

       // elastic_acceleration(dv, k, particles[tt-1], N);

        // Velocities / momenta
        x[2*nn] = - p->k * particles[tt-1][nn];

        // positions
        x[nn] += dv[nn];

        particles[tt][nn] = particles[tt-1][nn] + (dr + velocity_particles[nn]) * p->dt;      
    };
}
*/

int main(int argc, char *argv[]) {
    // char OUT_FILE [] = "plots.db";
    // char RUN_NAME[] = "elastic";
    std::cout << "starting... ";
   
    const Params * p = parse_cmd_line(argc, &*argv);
   
    cout << "variance : " << p->VARIANCE << endl;
   
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
        grid[0 ][ map_x_to_grid(particles[0][nn], p->GRID_X_POINTS) ]++;
    }
    
/*    for (size_t gg=0; gg < p->GRID_X_POINTS ; gg++){
        std::cout << grid[0][gg] << "  ";
    }
    cout << endl;
*/
  
    size_t tt = 0;
    print_to_db( db, p->RUN_NAME_STR.c_str(), grid[tt], tt, p->GRID_X_POINTS);
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

    float theta;
    float dr;
    std::normal_distribution<> normal_dist(0, p->VARIANCE);
   
    std::cout << "initialization completed." << std::endl;
    
    cout << grid[tt][0] << "\t" ;
    for (size_t tt=1; tt < p->T; tt++){
        for (size_t nn=0; nn < p->N ; nn++){
            // theta = std::generate_canonical<double, p->RANDOMNESS_BITS_U>(rnd2) ;

            dr = normal_dist(rnd1);

           // elastic_acceleration(dv, k, particles[tt-1], N);

            dv[nn] = - p->k * particles[tt-1][nn];

            velocity_particles[nn] += dv[nn] * p->dt;

            particles[tt][nn] = particles[tt-1][nn] + (dr + velocity_particles[nn]) * p->dt;
            
            int grid_point = map_x_to_grid(particles[tt][nn], p->GRID_X_POINTS);

            grid[tt][map_x_to_grid(particles[tt][nn], p->GRID_X_POINTS)]++; 
        };
       
        print_to_db( db, p->RUN_NAME_STR.c_str(), grid[tt], tt , p->GRID_X_POINTS);
        if (tt % 100 == 0){
            cout << '\r' <<  setfill(' ') << setw(5)  << tt;
        }
    }

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


