/**
 * @file mtcorr.h
 * @brief 
 * @version 0.1
 * @date 2022-12-26
 * @author Sota Miyamoto
 * @update Shoma Fujii
 * 
 * @copyright Copyright (c) 2022
 * 
 * @ref https://aip.scitation.org/doi/pdf/10.1063/1.3491098
 * 
 * usage:
 *   // declare //
 *   multipletau::correlator corr;
 *   // main loop //
 *   for i in 0..N
 *     corr(val)
 *   // get results //
 *   std::vector<double> t = corr.get_time_vec(); // must be mutiply dt
 *   std::vector<double> f = corr.get_corr_vec();
 * 
 * option:
 *   corr.save_state(string prefix, h5io::h5fp fp); 
 *       // save to prefix/**
 *   corr.load_state(string prefix, h5io::h5fp fp);
 *       // load from prefix/**
 */

#ifndef __MTCORR_H__
#define __MTCORR_H__

#include <iostream>
#include <stdexcept>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>

//#include"h5io.hpp"

namespace multipletau{
    class correlator
    {
        int S;
        int m;
        int p;
        Eigen::Vector3d *D;
        double *C;
        int    *N;
        Eigen::Vector3d *A;
        int    *M;

        public:
            correlator(int S_=40, int m_=2, int p_=16);
            ~correlator();
            void operator()(const Eigen::Vector3d &omega, const int i=0);

            std::vector<double> get_time_vec();
            std::vector<double> get_corr_vec();

            //void save_state(std::string prefix,h5io::h5fp& fp);
            //void load_state(std::string prefix, h5io::h5fp& fp);

        private:
            void alloc();
            void destroy();
            inline int ac2(const int x,const int y){
                return x*p+y;
            }
    };
} // namespace

// --------------------------------------------------------------------------------

namespace multipletau{
    correlator::correlator(int S_,int m_,int p_)
    {
        S=S_;
        m=m_;
        p=p_;
        alloc();
    }

    correlator::~correlator(){
        destroy();
    }

    void correlator::operator()(const Eigen::Vector3d &omega, const int i)
    {
        // (1) push omega to Dij
        for(int j = p-1; j > 0; j--){
            D[ac2(i,j)] = D[ac2(i,j-1)];
        }
        D[ac2(i,0)] = omega;

        // (2) update correlation array
        int begin,end;
        if(i==0){
            begin = 0; end = p;
        } else {
            begin = p/m; end = p;
        }
        for(int j = begin; j < end; j++)
        {
            C[ac2(i,j)] += D[ac2(i,0)].dot(D[ac2(i,j)]);
            N[ac2(i,j)] += 1;
        }

        // (3) omega is added to accumulator 
        A[i] += omega;
        M[i] += 1;

        // (4) if Mi==m the aceraged acc Ai/m sent to next
        if( M[i]==m ){
            if( i+1 == S+1 ){
                throw std::runtime_error("correlator buffer over flow!");
            }
            operator()(A[i]/(double)m, i+1);
            A[i] = Eigen::Vector3d::Zero();
            M[i] = 0;
        }
        return;
    }

    std::vector<double> correlator::get_time_vec(){
        std::vector<double> tk;
        double mi = 1.0;
        for(int i = 0; i < S+1; i++)
        {
            int begin,end;
            if(i==0){
                begin = 0; end = p;
            } else {
                begin = p/m; end = p;
            }
            for(int j = begin; j < end; j++){
                if (N[ac2(i,j)] == 0) return tk;
                tk.emplace_back(j*mi);
            }
            mi *= m;
        }
        return tk;
    }  

    std::vector<double> correlator::get_corr_vec(){
        std::vector<double> fk;
        for(int i = 0; i < S+1; i++)
        {
            int begin,end;
            if(i==0){
                begin = 0; end = p;
            } else {
                begin = p/m; end = p;
            }
            for(int j = begin; j < end; j++){
                if (N[ac2(i,j)] == 0) return fk;
                fk.emplace_back(C[ac2(i,j)]/N[ac2(i,j)]);
            }
        }
        return fk;
    }  

    //void correlator::save_state(std::string prefix,h5io::h5fp& fp){
    //  fp.create_dataset(prefix+"/S",{1},&S);    
    //  fp.create_dataset(prefix+"/m",{1},&m);    
    //  fp.create_dataset(prefix+"/p",{1},&p);
    //  fp.create_dataset(prefix+"/D",{S+1,p},D);    
    //  fp.create_dataset(prefix+"/C",{S+1,p},C);    
    //  fp.create_dataset(prefix+"/N",{S+1,p},N);    
    //  fp.create_dataset(prefix+"/A",{S+1  },A);    
    //  fp.create_dataset(prefix+"/M",{S+1  },M);    
    //  return;
    //}
    //
    //void correlator::load_state(std::string prefix, h5io::h5fp& fp){
    //  fp.read_dataset(prefix+"/S",&S);
    //  fp.read_dataset(prefix+"/m",&m);
    //  fp.read_dataset(prefix+"/p",&p);
    //  destroy();
    //  alloc();
    //  fp.read_dataset(prefix+"/D",D);
    //  fp.read_dataset(prefix+"/C",C);
    //  fp.read_dataset(prefix+"/N",N);
    //  fp.read_dataset(prefix+"/A",A);
    //  fp.read_dataset(prefix+"/M",M);
    //  return;
    //}

    void correlator::alloc(){
        D = new Eigen::Vector3d [p*(S+1)];
        C = new double [p*(S+1)];
        N = new    int [p*(S+1)];
        A = new Eigen::Vector3d [S+1];
        M = new    int [S+1];
        for(int i = 0; i < p*(S+1); i++){
            D[i] = Eigen::Vector3d::Zero();
            C[i] = 0;
            N[i] = 0;
        }
        for(int i = 0; i < (S+1); i++){
            A[i] = Eigen::Vector3d::Zero();
            M[i] = 0;
        }
    }

    void correlator::destroy(){
        delete[] D;
        delete[] C;
        delete[] N;
        delete[] A;
        delete[] M;
    }
}// namespace

#endif
