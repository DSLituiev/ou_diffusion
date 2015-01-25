#include <iostream>
#include <string.h> // memset
#include <stdio.h>
#include <stdlib.h>     /* abs */
#include <math.h>       /* round, floor, ceil, trunc */
#include <random>
#include <iomanip>      // std::setprecision
#include <sqlite3pp.h>

using namespace std;
size_t GRID_X_POINTS;

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
    gg = (abs(x_) < GRID_X_POINTS / 2) ? ( GRID_X_POINTS / 2 - x_) : (x_<0 ? 0: GRID_X_POINTS);
    return gg;
    }

void init_sql_table(sqlite3pp::database * db, const char run_name []){
        std::string sql_drop( "DROP TABLE IF EXISTS ");
        std::string name( run_name );
        sql_drop = sql_drop + run_name;
        db->execute(sql_drop.c_str());

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

void print_to_db(sqlite3pp::database * db, char * run_name, unsigned short int * grid_t, size_t & tt ){
    char buffer[GRID_X_POINTS*2];
    
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

    std::cout << "printed time point " << tt << " points: " <<  GRID_X_POINTS*2 <<" status :" << blob_status << std::endl;

    std::cout << sql << std::endl;
    for (size_t ii=0; ii < GRID_X_POINTS*2; ii++) {
        std::cout << (int) buffer[ii] << " ";
    }
/*
    std::cout << "time point: " << tt << std::endl;
    for (size_t gg=0; gg<GRID_X_POINTS ; gg++){
        std::cout << setfill(' ') << setw(3) << grid_t[gg] << "  ";
    }
*/
    std::cout << std::endl;
}

int main(void) {
    char OUT_FILE [] = "plots.db";
    char RUN_NAME[] = "simple";
    std::cout << "starting... ";
    int const N = 300;
    int const T = 20;
    int const RANDOMNESS_BITS_U = 10;
    GRID_X_POINTS = 1000-1;
    float VARIANCE = .2;
    
    float particles[T][N];
    for (size_t tt=0; tt<T; tt++){
        memset(particles[tt], 0, sizeof(float) * N);
    }

    unsigned short int grid[T][GRID_X_POINTS];
    for (size_t tt=0; tt<T; tt++){
        memset(grid[tt], 0, sizeof(unsigned short int) * GRID_X_POINTS);
    }

    for (size_t nn=0; nn<N ; nn++){
        grid[0 ][ map_x_to_grid(particles[0][nn], GRID_X_POINTS) ]++;
        std::cout << grid[0][map_x_to_grid(particles[0][nn], GRID_X_POINTS)] << "\t";
    }
    cout << endl;


    for (size_t gg=0; gg<GRID_X_POINTS ; gg++){
        std::cout << setfill(' ') << setw(3) << grid[0][gg] << "  ";
    }
    cout << endl;
    cout << "===========================================";
    cout << endl;

    // Good random seed, good engine
    auto rnd1 = std::mt19937(std::random_device{}());
    auto rnd2 = std::mt19937(std::random_device{}());

    float theta;
    float dr;
    std::normal_distribution<> normal_dist(0, VARIANCE);

    std::cout << "initialization completed." << std::endl;

    // DATABASE INITIALISATION
    sqlite3pp::database * db;
    db  = new sqlite3pp::database( OUT_FILE ) ;
    sqlite3pp::transaction * xct = new sqlite3pp::transaction(*db);
    init_sql_table(db, RUN_NAME);

    for (size_t tt=1; tt<T; tt++){
        for (size_t nn=0; nn<N ; nn++){
            theta = std::generate_canonical<double, RANDOMNESS_BITS_U>(rnd2) ;
            dr = normal_dist(rnd1);
            particles[tt][nn] = particles[tt-1][nn] + dr;
            int grid_point = map_x_to_grid(particles[tt][nn], GRID_X_POINTS);
//            cout << grid_point << "\t" ;
            grid[tt][map_x_to_grid(particles[tt][nn], GRID_X_POINTS)]++; 
        };
       
        print_to_db( db, RUN_NAME, grid[tt], tt );
 
        }

    xct->commit(); 
    std::cout << "finished successfully";
    std::cout << std::endl;
    return 0;
}


