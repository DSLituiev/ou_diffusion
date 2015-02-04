import numpy as np
import sqlite3
import struct
import matplotlib.pyplot as plt

def readByteIntoArray(lineList, N):
    # assert( len(lineList) == (1 + 2*rangel)*N ), 'wrong blob length: %u' % len(lineList)
    " split into N int numpy sub-arrays within list intList"
    it = struct.iter_unpack('=H', lineList)
    # intList = [[] for Null in range(N)]
    arr = np.empty(N, dtype = int);
    for i, temp_char in enumerate(it):
        # print('ind %u\t arr %u' % (i, i//(2*rangel+1)) )
        # intList[i].append(temp_char[0])
        arr[i] = temp_char[0]
    return arr


class SeriesSqlite():

    def __init__( self, db_name,  table_name, N):
        self.flag_read = False
        self.db_name = db_name
        self.table_name = table_name
        self.N = N
        self.track_num = 3
        
    def read_register(self):
        sql = 'SELECT * FROM register WHERE id = "%s"' % self.table_name;
        sql_params = 'SELECT grid_points, dt, grid_dx FROM register WHERE id = "%s"' % self.table_name;
        with sqlite3.connect(self.db_name) as  self.conn:
            curs = self.conn.cursor()
            curs.execute(sql)
            out = curs.fetchall()
            curs.execute(sql_params)
            out_params = curs.fetchall()
        self.N  = out_params[0][0]
        self.dt = out_params[0][1]
        self.dx = out_params[0][2]
        return out[0]
        
    def read(self, each = 0, from_ = 0, to_ = 0):
        self.read_register()
        
        conditions = []
        if each:
            conditions.append(' time %% %u == 0' %  each)
        if from_:
            conditions.append(' time > %u' % int(from_))
        if to_:
            conditions.append(' time < %u' % int(to_) )

        if not (each or from_ or to_):
            sql = 'SELECT * FROM %s' % self.table_name;
        else:
            sql = 'SELECT * FROM %s WHERE ' % self.table_name;
            sql += ' AND '.join(conditions)
            
        
        print(sql)
        
        with sqlite3.connect(self.db_name) as  self.conn:
            curs = self.conn.cursor()
            curs.execute(sql)
            out = curs.fetchall()
            Tp = len(out)
            
            out_arr = np.empty( (Tp, self.N), dtype = int )
            tracks =  np.empty( (Tp, self.track_num), dtype = float )
            t = np.empty( (Tp,1), dtype = float )
            
            for ti, line in enumerate(out):
                out_arr[ti, :] = readByteIntoArray(line[1], self.N)
                t[ti] = line[0]
                tracks[ti] = np.array(line[2:5])
            
            T = t[-1]
        self.flag_read = True
        return out_arr, Tp, t*self.dt , tracks, T
        
    def plot_data(self, *args, **kwargs):
        if not self.flag_read:
            out, Tp, t, tracks, T = ss.read(*args, **kwargs)        
        
        x = np.arange(-N//2,N//2) * self.dx
        fig = plt.figure()
        # fig, ax = plt.subplots()
        p = plt.pcolor(t[1:].T, x[1:-1], out[1:Tp,1:-1].T, vmin = 0, vmax = 100)
        fig.colorbar(p, shrink=0.5, aspect=8)
        plt.plot(t, tracks[:,0], 'kx-')
        # plt.plot(t, tracks[:,1], 'ko-')
        plt.plot(t, tracks[:,2], 'k.-')
        plt.title(table)
        plt.show()
        
        fig2 = plt.figure()
        plt.plot(t, tracks[:,0], 'kx-')
        plt.plot(t, tracks[:,1], 'ko-')
        plt.plot(t, tracks[:,2], 'k.-')
        plt.title(table)
        plt.show()
        return T
        
plt.close('all')
N = 400\

table =  "elastic006" # "simple" #
ss = SeriesSqlite("../bin/rw_simul-../plots.db", table, N)
info = ss.read_register()
print(info)

T = ss.plot_data(each = 50)
# print(out)
print( out[:, ss.N//2- 3:ss.N//2 + 3] )

# read again

out, Tp, t, tracks, _ = ss.read(each = 0, from_ = T*2/3 )

stat = np.mean(out[:, 1:-1], axis = 0)
x_clipped = x[1:-1]

fig = plt.figure()
plt.plot(x_clipped, stat)
plt.title(table)
plt.show()
 
mu = np.mean(stat*x_clipped) 
print('mean     : %f' % mu)
print('stdev : %f' %  np.sqrt(np.mean(stat*x_clipped*x_clipped - mu)) )

#fig = plt.figure()
#plt.plot(np.arange(0,Tp), 2*out[:,-1] / (out[:,-1]  + out[:,0]) - 1, 'rx')
#plt.plot(np.arange(0,Tp),
#         np.sum(out[:,1:-1] - stat, axis = 1)/ (N-2), 'bo')
#plt.show()


#fig = plt.figure()
#plt.plot(np.arange(0,Tp),
#         np.sum(out[:,1:-1] - stat, axis = 1)/ (N-2), 'bo')
#plt.show()

#fig = plt.figure()
#plt.plot(np.arange(0,Tp), (out[:,-1]  + out[:,0]) , 'rx')
#plt.show()

#fig = plt.figure()
#plt.plot(np.arange(0,Tp), out[:,-1] / (out[:,-1]  + out[:,0]) , 'rx')
#plt.show()

#fig = plt.figure()
#plt.plot( np.log10(x_clipped[N//2:]), np.log10(stat[N//2:]), '.-')
#plt.show()
#
#fig = plt.figure()
#plt.plot( np.log10(x_clipped[N//2:]),stat[N//2:], '.-')
#plt.show()

