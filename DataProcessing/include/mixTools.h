#ifndef MIXTOOLS
#define MIXTOOLS
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TMath.h"
#include "thrustTools.h"
#include "boostedEvtData.h"
#include "particleData.h"
#include "eventData.h"
#include <iostream>

//resets a mixed event's nParticle
void resetMixEvt(particleData *out){
  out->nParticle = 0;
}

void resetMixEvtBoosted(boostedEvtData *out){
  out->nParticle = 0;
}

//adds the particles from 'in' to tree 'out'
void appendMixEvt(particleData *out, particleData *in, TVector3 thrustAxis, TVector3 thrustAxis_ch, float nEventsFound){
  //setting nParticle correctly
  int oldnParticle = out->nParticle;
  out->nParticle += in->nParticle;
    
  //append extra particles to arrays
  for(int i = 0; i<in->nParticle; i++){
    out->px[oldnParticle+i] = in->px[i];
    out->py[oldnParticle+i] = in->py[i];
    out->pz[oldnParticle+i] = in->pz[i];
    out->pt[oldnParticle+i] = in->pt[i];
    out->pmag[oldnParticle+i] = in->pmag[i];
    out->rap[oldnParticle+i] = in->rap[i];
    out->eta[oldnParticle+i] = in->eta[i];
    out->theta[oldnParticle+i] = in->theta[i];
    out->phi[oldnParticle+i] = in->phi[i];
    out->mass[oldnParticle+i] = in->mass[i];
    out->charge[oldnParticle+i] = in->charge[i];
    out->pwflag[oldnParticle+i] = in->pwflag[i];
    out->pid[oldnParticle+i] = in->pid[i];
    out->d0[oldnParticle+i] = in->d0[i];
    out->z0[oldnParticle+i] = in->z0[i];
    out->ntpc[oldnParticle+i] = in->ntpc[i];
    out->nitc[oldnParticle+i] = in->nitc[i];
    out->nvdet[oldnParticle+i] = in->nvdet[i];
    out->vx[oldnParticle+i] = in->vx[i];
    out->vy[oldnParticle+i] = in->vy[i];
    out->vz[oldnParticle+i] = in->vz[i];
    
    TVector3 particle = TVector3(in->px[i],in->py[i],in->pz[i]);

    out->pt_wrtThr[oldnParticle+i] = ptFromThrust(thrustAxis, particle);
    out->eta_wrtThr[oldnParticle+i] = etaFromThrust(thrustAxis, particle);
    out->rap_wrtThr[oldnParticle+i] = rapFromThrust(thrustAxis, particle, in->mass[i]);
    out->theta_wrtThr[oldnParticle+i] = thetaFromThrust(thrustAxis, particle);
    out->phi_wrtThr[oldnParticle+i] = phiFromThrust(thrustAxis, particle);
    out->pt_wrtThrPerp[oldnParticle+i] = ptFromThrust(thrustAxis, particle, true);
    out->eta_wrtThrPerp[oldnParticle+i] = etaFromThrust(thrustAxis, particle, true);
    out->rap_wrtThrPerp[oldnParticle+i] = rapFromThrust(thrustAxis, particle, in->mass[i], true);
    out->theta_wrtThrPerp[oldnParticle+i] = thetaFromThrust(thrustAxis, particle, true);
    out->phi_wrtThr[oldnParticle+i] = phiFromThrust(thrustAxis, particle, true);
    out->pt_wrtChThr[oldnParticle+i] = ptFromThrust(thrustAxis_ch, particle);
    out->eta_wrtChThr[oldnParticle+i] = etaFromThrust(thrustAxis_ch, particle);
    out->rap_wrtChThr[oldnParticle+i] = rapFromThrust(thrustAxis_ch, particle, in->mass[i]);
    out->theta_wrtChThr[oldnParticle+i]= thetaFromThrust(thrustAxis_ch, particle);          
    out->phi_wrtChThr[oldnParticle+i] = phiFromThrust(thrustAxis_ch, particle);              
    out->pt_wrtChThrPerp[oldnParticle+i] = ptFromThrust(thrustAxis_ch, particle,true);
    out->eta_wrtChThrPerp[oldnParticle+i] = etaFromThrust(thrustAxis_ch, particle,true);
    out->rap_wrtChThrPerp[oldnParticle+i] = rapFromThrust(thrustAxis_ch, particle, in->mass[i],true);
    out->theta_wrtChThrPerp[oldnParticle+i]= thetaFromThrust(thrustAxis_ch, particle,true);          
    out->phi_wrtChThrPerp[oldnParticle+i] = phiFromThrust(thrustAxis_ch, particle,true);              
  }

  //account for fact that we are using multiple events to mix
  out->particleWeight = 1.0/nEventsFound;
}


//adds the particles from 'in' to tree 'out'
void appendMixEvtBoosted(boostedEvtData *out, particleData *in, TVector3 mixAxis, TVector3 mixBoost, float nEventsFound){
  //setting nParticle correctly
  int oldnParticle = out->nParticle;
  out->nParticle += in->nParticle;
    
  //append extra particles to arrays
  for(int i = 0; i<in->nParticle; i++){
    TLorentzVector particle = TLorentzVector(in->px[i],in->py[i],in->pz[i],in->mass[i]);
    particle.Boost(mixBoost);
    out->pt[oldnParticle+i] = ptFromThrust(mixAxis,particle.Vect());
    out->pmag[oldnParticle+i] = particle.Vect().Mag();
    out->eta[oldnParticle+i] = etaFromThrust(mixAxis,particle.Vect());
    out->theta[oldnParticle+i] = thetaFromThrust(mixAxis,particle.Vect());
    out->phi[oldnParticle+i] = phiFromThrust(mixAxis,particle.Vect());

    out->pt_Perp[oldnParticle+i] = ptFromThrust(mixAxis,particle.Vect(), true);
    out->pmag_Perp[oldnParticle+i] = particle.Vect().Mag();
    out->eta_Perp[oldnParticle+i] = etaFromThrust(mixAxis,particle.Vect(), true);
    out->theta_Perp[oldnParticle+i] = thetaFromThrust(mixAxis,particle.Vect(), true);
    out->phi_Perp[oldnParticle+i] = phiFromThrust(mixAxis,particle.Vect(), true);
  }

  //account for fact that we are using multiple events to mix
  out->particleWeight = 1.0/nEventsFound;
}

#endif

