! SPDX-FileCopyrightText: FLEXPART 1998-2019, see flexpart_license.txt
! SPDX-License-Identifier: GPL-3.0-or-later

function raerod (l,ust,z0)

  !*****************************************************************************
  !                                                                            *
  !     Calculation of the aerodynamical resistance ra from ground up to href  *
  !                                                                            *
  !     AUTHOR: Matthias Langer, modified by Andreas Stohl (6 August 1993)     *
  !                                                                            *
  !     Literature:                                                            *
  !     [1]  Hicks/Baldocchi/Meyers/Hosker/Matt (1987), A Preliminary          *
  !             Multiple Resistance Routine for Deriving Dry Deposition        *
  !             Velocities from Measured Quantities.                           *
  !             Water, Air and Soil Pollution 36 (1987), pp.311-330.           *
  !     [2]  Scire/Yamartino/Carmichael/Chang (1989),                          *
  !             CALGRID: A Mesoscale Photochemical Grid Model.                 *
  !             Vol II: User's Guide. (Report No.A049-1, June, 1989)           *
  !                                                                            *
  !     Variable list:                                                         *
  !     L     = Monin-Obukhov-length [m]                                       *
  !     ust   = friction velocity [m/sec]                                      *
  !     z0    = surface roughness length [m]                                   *
  !     href  = reference height [m], for which deposition velocity is         *
  !             calculated                                                     *
  !                                                                            *
  !     Constants:                                                             *
  !     karman    = von Karman-constant (~0.4)                                 *
  !     ramin = minimum resistence of ra (1 s/m)                               *
  !                                                                            *
  !     Subprograms and functions:                                             *
  !     function psih (z/L)                                                    *
  !                                                                            *
  !*****************************************************************************

  use par_mod

  implicit none

  real :: l,psih,raerod,ust,z0

  raerod=(alog(href/z0)-psih(href,l)+psih(z0,l))/(karman*ust)

end function raerod
