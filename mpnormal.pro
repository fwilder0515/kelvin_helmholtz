; ----- MPNORMAL ----- computes magnetopause normal for
;		        Fairfield magnetopause
;	Input:
;		Array of 3 elements containing the Position of S/C, 
;		in GSE and in Re
;	Output:
;		the Fairfield [1971] average magnetopause normal
;
;
;	Created by: Tai 95-09-28 (Adapted from Norbert Sckopke's xmpnormal)
;	Last modified: Tai 95-09-28


pro mpnormal, GSE, nvec

nvec= fltarr(3)


         PHI = ATAN(GSE(1),GSE(2))                                    
         RHO = SQRT (GSE(1)^2+GSE(2)^2)                               
         PSI = ACOS ((GSE(0)+25.3)/SQRT((GSE(0)+25.3)^2+8.02*RHO^2)) 

         nvec(0) = COS(PSI)
         nvec(1) = SIN(PSI) * SIN(PHI) 
         nvec(2) = SIN(PSI) * COS(PHI) 

;print,'normal=',nvec

;angle=xvangd(nvec,[1,0,0])

;print,'angle with sun-earth line',angle


RETURN
END
