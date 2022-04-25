#pragma once
#include "extra.h"
#include "energy.h"


template <class Tpsys>
PS::S32 UpdateRandM( Tpsys & pp,
		     PS::F64 time,
		     Extra*& cas,
		     PS::F64 & edisp,
		     PS::F64 mag) 
{
	    const PS::F64 cof = log(0.1);
	    const PS::F64 jor = 0.000284601;
            const PS::S32 n_loc = pp.getNumberOfParticleLocal();
	    PS::F64 edisp_loc = 0.0;
            PS::S32 cube = 0;
	    std::string s;
	    PS::F64 s_pha;
 
            for(PS::S32 i=0; i<n_loc; i++){
		s = "10000star";
		PS::S32 h = pp[i].ku;
		pp[i].pre_mass = pp[i].mass;
	        pp[i].pre_r_planet = pp[i].r_planet;	
		
    		while(1){
			  if( 10.0 <= pp[i].s_type ){
		              break;
			  }

			  cas[h].line = 0;
			  cas[h].setValue(s, h);

			  if( (cas[h].x1 <= time) && (time <= cas[h].x2) ){
			       pp[i].s_type = cas[h].pre_ty;
			       cas[h].interp1dim(time);
			       pp[i].r_planet = mag*cas[h].r;
			       pp[i].mass = cas[h].m;
			       break;
		           }
         		       cas[h].w = cas[h].w + 1;
	        }
		pp[i].delta_mass = pp[i].pre_mass - pp[i].mass;
             }
           PS::Comm::barrier();

            for(PS::S32 i=0; i<n_loc; i++){
		if(pp[i].delta_mass){
		   edisp_loc -= 0.5 * pp[i].delta_mass * pp[i].vel * pp[i].vel;
		   edisp_loc -= pp[i].delta_mass * pp[i].phi_s;
		   edisp_loc -= pp[i].delta_mass * pp[i].phi_d;
		   edisp_loc -= pp[i].delta_mass * pp[i].phi;
	        
		   cube = 1;
		}
#pragma omp parallel for reduction (-:edisp_loc)
		for(PS::S32 j=0; j<i; j++){
	            PS::F64 massi = pp[i].delta_mass;
	            PS::F64 massj = pp[j].delta_mass;
	            PS::F64vec posi = pp[i].pos;
		    PS::F64vec posj = pp[j].pos;
		    PS::F64vec dr = posi - posj;
		    PS::F64 eps2 = FP_t::eps2;

		    PS::F64 rinv = 1./sqrt(dr*dr + eps2);
		    edisp_loc -= massi * massj * rinv;
		}
            }
            PS::Comm::barrier();
            edisp += PS::Comm::getSum(edisp_loc);

    return cube;
}


template <class Tpsys>
void outputChange( Tpsys & pp,
		   PS::F64 time_sys,
		   std::ofstream & fout_chg,
		   PS::F64 de,
		   PS::F64 de_d_cum,
		   Energy e_init,
		   Energy e_now,
		   bool bSnap=true)
{
	if(PS::Comm::getRank() == 0 && bSnap){
           fout_chg << std::fixed << std::setprecision(8)
		    << time_sys << "\t"
		    << std::scientific << std::setprecision(15)
		    << de << "\t" << de_d_cum << "\t"
		    << e_init.etot << "\t" << e_now.etot << "\t" << e_now.edisp 
		    << std::endl;
	}
}

std::string setFname( std::string s,
		       PS::S32 i)
{
	std::string num = std::to_string(i+1);
	s = s + num + ".dat";
	return s;
}

PS::F64 leapEvolve( std::string s,  
		    PS::S32 k,
		    PS::S32 star_info)
{
   PS::F64 a;
   PS::F64 mom;
   PS::F64 phase;
   PS::F64 minfo;
   PS::F64 rinfo;
   PS::S32 f_counter = 1;
   PS::S32 s_counter = 1;
   PS::S32 t_counter = 1;

   std::ifstream ifs;
   ifs.open(setFname(s, k));
   if(ifs.fail()){
      std::cerr << "Can't file open\n";
      exit(1);
   }
   std::string sel;
   while(getline(ifs, sel)){
	 std::string sai;
	 std::istringstream stream(sel);
	 while(getline(stream, sai, ' ')){
	       a = stod(sai);
	       if(f_counter < 2){
		  mom = a;
		  f_counter++;
	       }
	       if(s_counter < 3){
		  phase = a;
		  s_counter++;
	       }
	       if(t_counter < 4){
		  minfo = a;
		  t_counter++;
	       }
	 }
	 f_counter = s_counter = t_counter = 1;
	 if(phase >= 10.0) break;
   }
   ifs.close();
   rinfo = a;

   if(star_info == 0){
      return rinfo;
   }else if(star_info == 1){
      return minfo;
   }else{
      return phase;
   }
}



template <class Tpsys>
PS::S32 ChangeAfterCol( Tpsys & pp,
		        PS::S32 n_col,
		        PS::F64 mag)
{
   const PS::S32 n_loc = pp.getNumberOfParticleLocal();
   PS::F64 edisp_loc = 0.0;

   std::string s = "10000star";

   PS::S32 col_type;
   PS::S32 tar_type;
   PS::S32 imp_type;

   PS::F64 pre_mass_tar = 0.;
   PS::F64 pre_mass_imp = 0.;
   PS::F64vec vel_tar;
   PS::F64vec vel_imp;
   PS::F64vec pos_tar;
   PS::F64vec pos_imp;

   PS::F64 imp_mass = 0.;
   PS::F64 tar_mass = 0.;

#pragma omp parallel for
   for(PS::S32 i=0; i<n_loc; i++){
       if( pp[i].tar_flag == 1 ){
	  tar_type = pp[i].s_type;
       }

       if( pp[i].imp_flag == 1 ){
	  imp_type = pp[i].s_type;
       }
   }

   PS::Comm::barrier(); 

//CollisionType
//nonGiant & nonGiant
   if( (tar_type<=1 || 7<=tar_type) && (imp_type<=1 || 7<=imp_type) ){
      col_type = 1;
      for(PS::S32 i=0; i<n_loc; i++){
          if(pp[i].imp_flag == 1){
            pp[i].imp_flag = 2;
            imp_mass = pp[i].mass;
	  }
	  if(pp[i].tar_flag == 1){
            pp[i].tar_flag = 2;
            tar_mass = pp[i].mass;
	  }
      }
   }
//nonGiant&Giant
   if( (2<=tar_type && tar_type<=6) && (imp_type<=1 || 7<=imp_type) ){
      col_type = 2;
      for(PS::S32 i=0; i<n_loc; i++){
         if(pp[i].imp_flag == 1){
            pp[i].imp_flag = 2;
            pp[i].tmp_imp_mass = pp[i].mass;
	    pp[i].r_planet = pp[i].pre_r_planet;
         }
	 if(pp[i].tar_flag == 1){
            pp[i].tar_flag = 2;
	    PS::S32 num = pp[i].ku;
            pp[i].r_planet = mag*leapEvolve(s, num, 0);
	    pp[i].tmp_tar_mass = pp[i].mass;
	    pp[i].mass = leapEvolve(s, num, 1);
	    pp[i].s_type = leapEvolve(s, num, 2);
	 }
      }
   }
//Giant&nonGiant
   if( (tar_type<=1 || 7<=tar_type) && (2<=imp_type && imp_type<=6) ){
      col_type = 3;
      for(PS::S32 i=0; i<n_loc; i++){
	  if(pp[i].imp_flag == 1){
             pp[i].imp_flag = 2;
             PS::S32 num = pp[i].ku;
	     pp[i].r_planet = mag*leapEvolve(s, num, 0);
	     pp[i].tmp_imp_mass = pp[i].mass;
	     pp[i].mass = leapEvolve(s, num, 1);
	     pp[i].s_type = leapEvolve(s, num, 2);
	  }
          if(pp[i].tar_flag ==1){
             pp[i].tar_flag = 2;
             pp[i].tmp_tar_mass = pp[i].mass;
	     pp[i].r_planet = pp[i].pre_r_planet;
          }
      }
   }
//Giant&Giant
   if( (2<=tar_type && tar_type<=6) && (2<=imp_type && imp_type<=6) ){
      col_type = 4;
      for(PS::S32 i=0; i<n_loc; i++){
	  if(pp[i].imp_flag == 1){
             pp[i].imp_flag = 2;
	     PS::S32 num = pp[i].ku;
	     pp[i].r_planet = mag*leapEvolve(s, num, 0);
	     pp[i].tmp_imp_mass = pp[i].mass;
	     pp[i].mass = leapEvolve(s, num, 1);
	     pp[i].s_type = leapEvolve(s, num, 2);
	  }
          if(pp[i].tar_flag == 1){
             pp[i].tar_flag = 2;
	     PS::S32 num = pp[i].ku;
	     pp[i].r_planet = mag*leapEvolve(s, num, 0);
	     pp[i].tmp_tar_mass = pp[i].mass;
	     pp[i].mass = leapEvolve(s, num, 1);
	     pp[i].s_type = leapEvolve(s, num, 2);
	  }
      }
   }

   std::cout << "col_type: " << col_type << std::endl;
   return col_type;
}


template <class Tpsys>
void colTypeTtoF( Tpsys & pp,
		  PS::F64 & edisp)
{
      const PS::S32 n_loc = pp.getNumberOfParticleLocal();
      PS::F64 edisp_loc = 0.;
      PS::F64 imp_change;
      PS::F64 tar_change;
      PS::F64vec pos_tar;
      PS::F64vec pos_imp;

      for( PS::S32 i=0; i<n_loc; i++ ){
          PS::F64 tar_tot = 0.;
          PS::F64 imp_tot = 0.;
          if( pp[i].tar_flag==2 ){
             pp[i].delta_mass_tar = pp[i].tmp_tar_mass - pp[i].mass;
             edisp_loc -= 0.5 * pp[i].delta_mass_tar * pp[i].vel * pp[i].vel;
             edisp_loc -= pp[i].delta_mass_tar * pp[i].phi_s;
             edisp_loc -= pp[i].delta_mass_tar * pp[i].phi_d;
             edisp_loc -= pp[i].delta_mass_tar * pp[i].phi;
	     tar_change = pp[i].delta_mass_tar;
	     pos_tar = pp[i].pos;
             pp[i].tar_flag = 0;
          }
          if( pp[i].imp_flag==2 ){
	     pp[i].delta_mass_imp = pp[i].tmp_imp_mass - pp[i].mass;
             edisp_loc -= 0.5 * pp[i].delta_mass_imp * pp[i].vel * pp[i].vel;
             edisp_loc -= pp[i].delta_mass_imp * pp[i].phi_s;
             edisp_loc -= pp[i].delta_mass_imp * pp[i].phi_d;
             edisp_loc -= pp[i].delta_mass_imp * pp[i].phi;
	     imp_change = pp[i].delta_mass_imp;
	     pos_imp = pp[i].pos;
             pp[i].imp_flag = 0;
          }
      }
        PS::Comm::barrier();

	PS::F64 eps2 = FP_t::eps2;
	PS::F64vec dr = pos_tar - pos_imp; 
	PS::F64 rinv = 1.0/sqrt(dr*dr+eps2);
	edisp_loc -= tar_change * imp_change * rinv;

        edisp += PS::Comm::getSum(edisp_loc);

}
  
/*
   if(col_type == 1){
      for(PS::S32 i=0; i<n_loc; i++){
	  if(pp[i].tar_flag == 1){
             remove_tar[n_remove_tar] = i;
	     n_remove_tar++; 
             edisp_loc -= 0.5 * pp[i].mass * pp[i].vel * pp[i].vel;
	     edisp_loc -= pp[i].mass * pp[i].phi_s;
	     edisp_loc -= pp[i].mass * pp[i].phi_d;
	     edisp_loc -= pp[i].mass * pp[i].phi;
	     tar_change = pp[i].mass;
	     pos_tar = pp[i].pos;
	  }
	  if(pp[i].imp_flag == 1){
             remove_imp[n_remove_imp] = i;
	     n_remove_imp++;
      	     edisp_loc -= 0.5 * pp[i].mass * pp[i].vel * pp[i].vel;
	     edisp_loc -= pp[i].mass * pp[i].phi_s;
	     edisp_loc -= pp[i].mass * pp[i].phi_d;;
	     edisp_loc -= pp[i].mass * pp[i].phi;
	     imp_change = pp[i].mass;
	     pos_imp = pp[i].pos;
	  }
      }

	PS::F64 eps2 = FP_t::eps2;
	PS::F64vec dr = pos_tar - pos_imp; 
	PS::F64 rinv = 1.0/sqrt(dr*dr+eps2);
	edisp_loc += - tar_change * imp_change * rinv;

	
	edisp += PS::Comm::getSum(edisp_loc);
     }


      //fout_colrem << "Collision & Remove Particle" << "\n" 
      fout_colrem << std::scientific << std::setprecision(15) << tar_change <<"\t"
		  << imp_change
		  << std::endl;
*/
template <class Tpsys>
void colTypeO( Tpsys & pp,
               PS::F64 & edisp,
               PS::S32 n_col)
{
     const PS::S32 n_loc = pp.getNumberOfParticleLocal();
     PS::S32 n_colrem = 0;
     PS::S32 s = n_col * 2;
     PS::S32 t = 2;
     PS::S32 * colrem = new PS::S32[t];
     PS::F64 edisp_loc = 0.;
     PS::F64 t_mass = 0.;
     PS::F64 i_mass = 0.;
     PS::F64vec pot;
     PS::F64vec poi;
     //PS::F64 tar_tot = 0.;
     //PS::F64 imp_tot = 0.;

     for(PS::S32 i=0; i<n_loc; i++){
         //PS::F64 tar_tot = 0.;
         //PS::F64 imp_tot = 0.;
         if( pp[i].tar_flag == 2 ){
             edisp_loc -= 0.5 * pp[i].mass * pp[i].vel * pp[i].vel;
             edisp_loc -= pp[i].mass * pp[i].phi_s;
             edisp_loc -= pp[i].mass * pp[i].phi_d;
             edisp_loc -= pp[i].mass * pp[i].phi;
             pot = pp[i].pos;
             t_mass = pp[i].mass;
             //{
                    colrem[0] = i;
                    n_colrem++;
             //}
         }
         if( pp[i].imp_flag == 2 ){
             edisp_loc -= 0.5 * pp[i].mass * pp[i].vel * pp[i].vel;
             edisp_loc -= pp[i].mass * pp[i].phi_s;
             edisp_loc -= pp[i].mass * pp[i].phi_d;
             edisp_loc -= pp[i].mass * pp[i].phi;
            //imp_tot = i_kin + i_phi_s + i_phi_d + i_phi;
             poi = pp[i].pos;
             i_mass = pp[i].mass;
             //{
                     colrem[1] = i;
                     n_colrem++;
             //}
         }
         //edisp_loc -= tar_tot + imp_tot;
         if( n_colrem == 2 ){
             PS::F64 eps2 = FP_t::eps2;
             PS::F64vec dr = pot - poi;
             PS::F64 rinv = 1./sqrt(dr*dr + eps2);
             edisp_loc += - t_mass * i_mass * rinv;
         }
     }
/*
     PS::F64 eps2 = FP_t::eps2;
     PS::F64vec dr = pot - poi;
     PS::F64 rinv = 1./sqrt(dr*dr + eps2);
     edisp_loc -= t_mass * i_mass * rinv;
*/
     edisp += PS::Comm::getSum(edisp_loc);

     if( n_colrem ){
             pp.removeParticle(colrem, n_colrem);
     }
     delete [] colrem;

#pragma omp prallel for
     for(PS::S32 i=0; i<n_loc; i++){
          pp[i].imp_flag = 0;
          pp[i].tar_flag = 0;
     }

}

template <class Tpsys>
void colTypeO( Tpsys & pp,
	       PS::F64 & edisp)
{
     const PS::S32 n_loc = pp.getNumberOfParticleLocal();
     const PS::S32 n_proc = PS::Comm::getNumberOfProc();
     PS::F64 edisp_loc = 0.;

     static std::vector<PS::S32> n_colrem_list;
     static std::vector<PS::S32> n_colrem_adr;
     static std::vector<FP_t> colrem_list_loc;
     static std::vector<FP_t> colrem_list_glb;
     n_colrem_list.resize(n_proc);
     n_colrem_adr.resize(n_proc);

     static std::vector<PS::S32> colrem_list;
     colrem_list.clear();

#pragma omp parallel for
       for(PS::S32 i=0; i<n_loc; i++){
	   if(pp[i].tar_flag==2 || pp[i].imp_flag==2){
#pragma omp critical		   
              {
	               colrem_list.push_back(i);
                       pp[i].tar_flag = 0;
                       pp[i].imp_flag = 0;
	      }
	   }
     }

     PS::S32 n_colrem_loc = colrem_list.size();
     PS::S32 n_colrem_glb = PS::Comm::getSum(n_colrem_loc);
    
     if( n_colrem_glb ){
	     
         if(PS::Comm::getRank() == 0){
                 colrem_list_glb.resize(n_colrem_glb);
         }
         colrem_list_loc.resize(n_colrem_loc);

#ifdef PARTICLE_SIMULATOR_PARALLEL
	 MPI_Gather(&n_colrem_loc, 1, PS::GetDataType(n_colrem_loc),
	            &n_colrem_list[0], 1, PS::GetDataType(n_colrem_list[0], 0, MPI_COMM_WORLD);
#else
         n_colrem_list[0] = n_colrem_loc;
#endif
	 if( PS::Comm::getRank() == 0 ){
	    PS::S32 tmp_colrem = 0;
	    for ( PS::S32 i=0; i<n_proc; i++){
	        n_colrem_adr[i] = tmp_colrem;
		tmp_colrem += n_colrem_list[i];
	    }
	    assert ( n_colrem_glb == tmp_colrem );
	 }

         for(PS::S32 i=0; i<n_colrem_loc; i++){
	         colrem_list_loc[i] = pp[colrem_list.at(i)];
         }

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
	 MPI_Gatherv(&colrem_list_loc[0], n_colrem_loc,                        PS::GetDataType(colrem_list_loc[0]),
		     &colrem_list_glb[0], &n_colrem_list[0], &n_colrem_adr[0], PS::GetDataType(colrem_list_glb[0]), 0, MPI_COMM_WORLD);
#else
         for(PS::S32 i=0; i<n_colrem_loc; i++) colrem_list_glb[i] = colrem_list_loc[i];
#endif

         if( PS::Comm::getRank() == 0 ){
            for(PS::S32 i=0; i<n_colrem_glb; i++){

	        PS::F64 massi = colrem_list_glb[i].mass;
		PS::F64vec veli = colrem_list_glb[i].vel;
	        edisp_loc -= 0.5 * massi * veli * veli;
	        edisp_loc -= massi * colrem_list_glb[i].phi_s;
	        edisp_loc -= massi * colrem_list_glb[i].phi_d;
	        edisp_loc -= massi * colrem_list_glb[i].phi;
	  
	        for(PS::S32 j=0; j<i; j++){
	           if(colrem_list_glb[j].id != colrem_list_glb[i].id){
		      PS::F64 massj = colrem_list_glb[j].mass;
		      PS::F64vec posi = colrem_list_glb[i].pos;
		      PS::F64vec posj = colrem_list_glb[j].pos;
		      PS::F64 eps2 = FP_t::eps2;

		      PS::F64vec dr = posi - posj;
		      PS::F64    rinv = 1./sqrt(dr*dr + eps2);

		      edisp_loc += - massi * massj * rinv;
	           }  	    
                }   
            }
         }
#ifdef INDIRECT_TERM
	 e_ind_before = calcIndirectEnergy(pp);
#endif
      }

      if (n_colrem_loc) pp.removeParticle(&colrem_list[0], n_colrem_loc);

      edisp += PS::Comm::getSum(edisp_loc);

/*
    for(PS::S32 i=0; i<n_loc; i++){
	pp[i].tar_flag = 0;
	pp[i].imp_flag = 0;
    } 
i*/

}






