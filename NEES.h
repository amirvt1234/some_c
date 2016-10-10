//    An input file example to run a case based on GPUSPH code

#include "Problem.h"
#include "Point.h"
#include "Cube.h"
#include "Rect.h"
#include "Cylinder.h"
#include "Vector.h"
#include "Cone.h"

class NEES: public Problem {
	private:
		int			icyl, icone, wmakertype;
		Cube		experiment_box;
		Rect		experiment_box1, experiment_box2;
		int			i_use_bottom_plane;
		PointVect	parts;
		PointVect	boundary_parts;
		PointVect	piston_parts;
		PointVect	gate_parts;

		double 		lx, ly, lz;	// Dimension of the computational domain
		double		h_length;	// Horizontal part of the experimental domain
		double		slope_length;	// Length of the inclined plane
		double		height;		// Still water (with z origin on the horizontal part)
		double		beta;		// Angle of the inclined plane
       	        double		H;		// still water level
		double		Hbox;	// height of experiment box

		// Moving boundary data
		double		m_S, m_Hoh, m_tau;

	public:
		NEES(const GlobalData *);
		~NEES(void);
		int fill_parts(void);
		uint fill_planes(void);
		void copy_planes(float4*, float*);

		void copy_to_array(BufferList &);
		MbCallBack& mb_callback(const float, const float, const int);

		void release_memory(void);
};

