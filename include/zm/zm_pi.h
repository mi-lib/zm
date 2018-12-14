/* ZM - Z's Mathematics Toolbox
 * Copyright (C) 1998 Tomomichi Sugihara (Zhidao)
 */
/*! \file zm_pi.h
 * \brief circle ratio.
 * \author Zhidao
 */

#ifndef __ZM_PI_H__
#define __ZM_PI_H__

/*! \brief circle ratio (pi). */
#define zPI   3.14159265358979323846
/*! \brief circle ratio times 2 (2 pi). */
#define zPIx2 6.28318530717958647692
/*! \brief circle ratio divided by 2 (pi / 2). */
#define zPI_2 1.57079632679489661923

/*! \brief radian per degree (pi / 180). */
#define zRAD_PER_DEG  0.01745329251994329547
/*! \brief degree per radian (180 / pi). */
#define zDEG_PER_RAD 57.29577951308232286465

/*! \brief convert a value \a x in degree to a radian value. */
#define zDeg2Rad(x) ( (x)*zRAD_PER_DEG )
/*! \brief convert a value \a x in radian to a degree value. */
#define zRad2Deg(x) ( (x)*zDEG_PER_RAD )

#endif /* __ZM_PI_H__ */
