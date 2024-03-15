   !---------------------------------------------------------------------------------------------

  subroutine Dike_Emplace(d_erupt_xy,T,T_melt,f,df,C,vm,nx,ny,nz,nb,zFS,sig_y,DikeN,EDikeN,RDikeN,sigma_dike,OverP,NewDike) 

   ! This module evaluates the overpressure required to form a dike and flags cells
   ! that meet or exceed the lithospheric yield stress for diking and eruption.
   ! It then calculates the penetration depth of the dike and feeds that into df_erupt_intrude.
     
     use grid,only: gl,dzg,D_dimensional
     use fundconst, only: pi
     use comppars,only: nCompFields
     use refstat, only: g_dimensional,grav0,Cp_dimensional,dens_dimensional,tcond_dimensional,fCp
     use meltstuff, only: ddens_sol_liq_dimensional,eta_magma_dimensional,K0_dimensional,latent_dike,erupt_T_lith,DefDikeDT,DTDike,CpDike,Set_vDike,vDike,DikeQDist,DikeQSeed,DikeQCent 
     use timestep, only: total_time
     use io, only: fname
     implicit none
     integer,intent(in):: nx,ny,nz,nb
     real,intent(in):: f(-1:nx,-1:ny,-1:nz,nb,nMeltFields), df(-1:nx,-1:ny,-1:nz,nb,nMeltFields)
     real,intent(in):: d_erupt_xy(-1:nx,-1:ny,nb),zFS(-1:nx,-1:ny,nb),T_Melt(-1:nx,-1:ny,nb),vm(-1:nx,-1:ny,-1:nz,nb)
     real,intent(in):: sig_y,T(-1:nx,-1:ny,-1:nz,nb), C(-1:nx,-1:ny,-1:nz,nb,nCompFields)
     real,intent(out):: sigma_dike(-1:nx,-1:ny,nb),OverP(-1:nx,-1:ny,-1:nz,nb)   
     real,intent(inout):: DikeN(-1:nx,-1:ny,-1:nz,nb),EDikeN(-1:nx,-1:ny,-1:nz,nb),RDikeN(-1:nx,-1:ny,-1:nz,nb)
     logical,intent(out):: NewDike(-1:nx,-1:ny,nb)
     real,allocatable:: h(:,:,:),hm(:,:,:),Tz(:),QrandArray(:,:,:)!,fze_new(:,:,:,:),fze_tot_new(:,:,:),Kze(:,:,:)
     !real,allocatable:: CpDike(:,:,:),T_erupt(:,:,:),C_erupt(:,:,:,:),f_erupt(:,:,:,:)
     real:: f_new(nMeltFields),f_tot_new,depth,DikeV,DikeQ,DikeW,DikeKappa,Diketc,DikeD,DikeDT,LambdaLHS,LambdaRHS,Lambda2,LambdaDike,TLAvg,DikeProp,fze_new(nMeltFields),fze_tot_new
     real:: SigOne,SigTwo,SigThree,SigFour,SigFive,SigSix,SigSev,Qrand, dDTdz, dDTdt, tDT, DTn, Diketcn, LambdaDiken, LambdaLHSn, DikeK
     integer,allocatable,dimension(:,:,:):: z_erupt
     integer ix,iy,iz,ib,ze,Lambda1, DTnum, n
     character(len=100) :: filename
 
 
     allocate (z_erupt(-1:nx,-1:ny,nb),h(-1:nx,-1:ny,nb),hm(-1:nx,-1:ny,nb),Tz(-1:nz),QrandArray(-1:nx,-1:ny,nb))!,fze_new(-1:nx,-1:ny,nb,1:nMeltFields),fze_tot_new(-1:nx,-1:ny,nb),Kze(-1:nx,-1:ny,nb))
     !allocate(CpDike(-1:nx,-1:ny,nb),T_erupt(-1:nx,-1:ny,nb),C_erupt(-1:nx,-1:ny,nb,nCompFields),f_erupt(-1:nx,-1:ny,nb,nMeltFields))
     h(:,:,:) = 0.0
     hm(:,:,:) = 0.0
     NewDike(:,:,:) = .false.
     OverP(:,:,:,:) = 0.0
     if (DikeQDist) then
         if (DikeQSeed) then
            call random_seed()
         end if
         call random_number(QrandArray)
     end if

     filename = trim(fname) // 'DikeStats.dat'
     open(66,file=filename,access='append')
     
     
 
    !Define z node id at base of lithosphere for all x, y
     do concurrent(ib=1:nb,iy=0:ny-1)
        do ix = gl%yy_ixmin(iy),gl%yy_ixmax(iy)
           do iz=0,nz-1
              depth = zFS(ix,iy,ib) - gl%z(1,iz)
              if(depth==d_erupt_xy(ix,iy,ib)) then
                 z_erupt(ix,iy,ib)=iz
                 !T_erupt(ix,iy,ib)=T(ix,iy,iz,ib)
                 !C_erupt(ix,iy,ib,:)=C(ix,iy,iz,ib,:)
                 !f_erupt(ix,iy,ib,:)=f(ix,iy,iz,ib,:) !+ df(ix,iy,iz,ib,:)
              end if
           end do
        end do
     end do

     
     !For now to keep things simple some dike physical parameters (e.g. specific heat, latent heat) are defined separately from wider StagYY scaling.
     !This can be fixed when merging with the official StagYY build.
     !CpDike(:,:,:)=fCp(z_erupt(:,:,:),1,T_erupt(:,:,:),T_erupt(:,:,:),T_erupt(:,:,:),C_erupt(:,:,:,:),f_erupt(:,:,:,:))
     !fze_new(:,:,:,:) = f(:,:,z_erupt(:,:,:),:,:) + df(:,:,z_erupt(:,:,:),:,:)
     !fze_tot_new(:,:,:) = sum(fze_new(:,:,:,1:nMeltFields))
     !Kze(:,:,:) = permeability_over_eta_m(fze_tot_new(:,:,:))
 
     do concurrent(ib=1:nb,iy=0:ny-1)
        do ix = gl%yy_ixmin(iy),gl%yy_ixmax(iy)
           ze = z_erupt(ix,iy,ib)
           iz = ze
           fze_new (:) = f(ix,iy,ze-1,ib,:) + df(ix,iy,ze-1,ib,:)
           fze_tot_new = sum(fze_new(1:nMeltFields))
           do while (iz>0) !determine height of magma in column
              iz = iz-1
              f_new(:)  = f(ix,iy,iz,ib,:) + df(ix,iy,iz,ib,:)
              f_tot_new = sum(f_new(1:nMeltFields))
              h(ix,iy,ib) = h(ix,iy,ib) + (gl%z(0,iz)-gl%z(0,iz-1))*f_tot_new
              hm(ix,iy,ib) = hm(ix,iy,ib) + (gl%z(0,iz)-gl%z(0,iz-1))
              if (f_tot_new == 0) exit
              cycle
           end do
             !determine buoyancy stress of magma, flag new dike if > sigy
           sigma_dike(ix,iy,ib) = ddens_sol_liq_dimensional * h(ix,iy,ib) * (g_dimensional * grav0)
           OverP(ix,iy,ze,ib) = sigma_dike(ix,iy,ib)  
           if (sigma_dike(ix,iy,ib) >= sig_y) then 
              DikeN(ix,iy,ze,ib) = DikeN(ix,iy,ze,ib) + 1. !Increment dike counter
              DikeK = permeability_over_eta_m(fze_tot_new)
              DikeV = (DikeK * (ddens_sol_liq_dimensional*grav0))*((1.-fze_tot_new)/fze_tot_new)
              !DikeV = vm(ix,iy,ze-1,ib)
              if (set_vDike) then
                 DikeV = vDike * DikeV
              !elseif (DikeV == 0) then
                 !DikeV = vm(ix,iy,(ze-1),ib)
              end if
              !Compute relevant dike properties to determine cooling time and lithospheric penetration
              if (DikeQDist) then   ! Assigns dike a pseduo-random normally distributed Qd value with a mean Qd of 0.1 m^2/s to simulate natural variation of dike magnitude
                  SigOne = 0.001
                  SigTwo = 0.022
                  SigThree = 0.158
                  SigFour = 0.5
                  SigFive = 0.842
                  SigSix = 0.978
                  SigSev = 0.999
                  Qrand = QrandArray(ix,iy,ib)
                  if (Qrand<=SigOne) then
                     DikeQ = (Qrand/SigOne) * DikeQCent * 1.e-3
                  else if (Qrand<=SigTwo) then
                     DikeQ = ((Qrand-SigOne)/(SigTwo-SigOne)) * DikeQCent * 1.e-2
                  else if (Qrand<=SigThree) then
                     DikeQ = ((Qrand-SigTwo)/(SigThree-SigTwo)) * DikeQCent * 1.e-1
                  else if (Qrand<=SigFour) then
                     DikeQ = ((Qrand-SigThree)/(SigFour-SigThree)) * DikeQCent
                  else if (Qrand<=SigFive) then
                     DikeQ = ((Qrand-SigFour)/(SigFive-SigFour)) * DikeQCent * 1.e1
                  else if (Qrand<=SigSix) then
                     DikeQ = ((Qrand-SigFive)/(SigSix-SigFive)) * DikeQCent * 1.e2
                  else if (Qrand<=SigSev) then
                     DikeQ = ((Qrand-SigSix)/(SigSev-SigSix)) * DikeQCent * 1.e3
                  else if (Qrand<=1.0) then
                     DikeQ = ((Qrand-SigSev)/(1.0-SigSev)) * DikeQCent * 1.e4
                  end if
               else
                  DikeQ = DikeV * hm(ix,iy,ib)
               end if
              DikeW = ((3.*eta_magma_dimensional*DikeQ)/(2.*ddens_sol_liq_dimensional*grav0))**(0.33333)
              DikeProp =((DikeQ**2*ddens_sol_liq_dimensional*grav0)/(12*eta_magma_dimensional))**(0.333333)
              DikeKappa = 3.0 / (3300. * CpDike) !placeholder for tcond_dimensional / (dens_dimensional*CpDike) !(ix,iy,ib))
              if (DefDikeDT) then
                DikeDT = DTDike
                LambdaLHS = (latent_dike * sqrt(pi))/(DikeDT*CpDike) !(ix,iy,ib))
                do Lambda1 = 0, 13000, 1
                   Lambda2 = Lambda1 * 0.0001
                   LambdaRHS = exp(-Lambda2**2)/(Lambda2*(1+erf(Lambda2)))
                   if (abs(LambdaLHS - LambdaRHS)<=0.001) then
                       LambdaDike = Lambda2
                       exit
                   end if
                end do
                Diketc = (DikeW**2)/(4*DikeKappa*LambdaDike)
              else
                !Tz = T(ix,iy,:,ib)
                !Tz(-1:(ze-1)) = 0.
                !TLAvg = (sum(Tz))/(nz-ze)
                !DikeDT = T_melt(ix,iy,ib) - TLAvg
                dDTdz = ((erupt_T_lith - 273.) / d_erupt_xy(ix,iy,ib)) * 1000.
                DTnum = int(d_erupt_xy(ix,iy,ib) / 1000.)
                DTiterate: do n = 1, DTnum, 1
                  DTn = n * dDtDz
                  dDTdt = (dDTdz / 1000.) * DikeProp
                  tDT = DTn / dDTdt
                  LambdaLHSn = (latent_dike * sqrt(pi))/(DTn*CpDike) !(ix,iy,ib))
                  Cooling: do Lambda1 = 0, 130000, 1
                     Lambda2 = Lambda1 * 0.00001
                     LambdaRHS = exp(-Lambda2**2)/(Lambda2*(1+erf(Lambda2)))
                     if ((LambdaLHSn >= 10.0).and.(abs(LambdaLHSn - LambdaRHS)<=1.0)) then
                         LambdaDiken = Lambda2
                         exit Cooling
                     else if ((LambdaLHSn>=1.0).and.(abs(LambdaLHSn - LambdaRHS)<=0.1)) then
                         LambdaDiken = Lambda2
                         exit Cooling
                     else if ((LambdaLHSn<1.0).and.(abs(LambdaLHSn - LambdaRHS)<=0.001)) then
                        LambdaDiken = Lambda2
                        exit Cooling
                     end if
                  end do Cooling
                  Diketcn = (DikeW**2)/(4*DikeKappa*LambdaDiken**2)
                  if (tDT > Diketcn) then
                     Diketc = Diketcn
                     LambdaDike = LambdaDiken
                     DikeDT = DTn
                     LambdaLHS = LambdaLHSn
                     exit DTiterate
                  else if (n==DTnum) then
                     Diketc = Diketcn
                     LambdaDike = LambdaDiken
                     DikeDT = DTn
                     LambdaLHS = LambdaLHSn
                  end if
                end do DTiterate

                     

              end if

              
              
              DikeD = Diketc * DikeProp
              if (DikeD>=D_erupt_xy(ix,iy,ib)) then
                 NewDike(ix,iy,ib) = .true.
                 EDikeN(ix,iy,ze,ib) = EDikeN(ix,iy,ze,ib) + 1. !Increment eruption counter
              end if
              RDikeN(ix,iy,ze,ib) = EDikeN(ix,iy,ze,ib) / DikeN(ix,iy,ze,ib)
              print*,'Dike triggered at coord x = ',ix,', y = ',iy,', z = ',ze
              print*,'DikeV       = ',DikeV,', DikeQ      = ',DikeQ,', DikeW   = ',DikeW
              print*,'etam        = ',eta_magma_dimensional,', ddens      = ',ddens_sol_liq_dimensional,', g       = ',grav0
              print*,'T_melt      = ',T_melt(ix,iy,ib),', TLavg      = ',TLavg,', DikeDT  = ',DikeDT
              print*,'latent_dike = ',latent_dike,', cp         = ',CpDike !(ix,iy,ib)
              print*,'LambdaLHS   = ',LambdaLHS,', LambdaDike = ',LambdaDike,', DikeProp= ',DikeProp
              print*,'Diketc      = ',Diketc,', DikeD      = ',DikeD,', d_erupt = ',d_erupt_xy(ix,iy,ib)
              print*,'DikeKappa   = ',DikeKappa,', tcond      = ',tcond_dimensional,', density = ',dens_dimensional
              print*,'fze_tot_new = ',fze_tot_new,', DikeK      = ',DikeK 
              write(66,'(I5,I5,I5,8es15.7,8es15.7,8es15.7,8es15.7,8es15.7,8es15.7,8es15.7,8es15.7)') ix,iy,ze,total_time,DikeV,DikeQ,DikeW,DikeProp,Diketc,DikeD,d_erupt_xy(ix,iy,ib)
           end if
        end do
     end do
    close(66)  
   end subroutine Dike_Emplace 

   !-------------------------------------------------------------------------------------------

  pure function d_erupt_adjust(f,df,nx,ny,nz,nb,xe,ye,ze,be,zFS) result (depth)
     
     use grid,only: gl,dzg,D_dimensional
     implicit none
     real depth
     integer,intent(in):: nx,ny,nz,nb,xe,ye,ze,be
     real,intent(in):: f(-1:nx,-1:ny,-1:nz,nb,nMeltFields), df(-1:nx,-1:ny,-1:nz,nb,nMeltFields)
     real,intent(in)::zFS(-1:nx,-1:ny,nb)
     real:: f_new(nMeltfields),f_tot_new, d
     integer iz


     do iz=nz-1,ze,-1
        d = zFS(xe,ye,be) - gl%z(1,ze)
        f_new = f(xe,ye,iz,be,:) + df(xe,ye,iz,be,:)
        f_tot_new = sum(f_new(1:nMeltFields))
        if (f_tot_new > 0) exit
     end do

     depth = d

  end function d_erupt_adjust
