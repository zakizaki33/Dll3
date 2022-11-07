#ifndef CXYLIST_H_INCLUDED
#define CXYLIST_H_INCLUDED

#include "list.h"
#include "cFitting.h"

class cXYData {
public:
	double x;
	double y;
	cXYData();
	cXYData(double x,double y=1);
	bool operator>(const cXYData& a) const;
	bool operator<(const cXYData& a) const;
	friend std::istream& operator>>(std::istream& from, cXYData& a);
	friend std::ostream& operator<<(std::ostream& to, const cXYData& a);
};

class cXYList : public list<cXYData>
{
public:
	void AddData(double x,double y);
	double x(int i);
	double y(double x);
	double dydx(double x);
	cXYList Subset(double x_start,double x_end);
	void Apply(cXYList &x);
	double Integral(double x_start,double x_end);
	double Average(double x_start,double x_end);
	double Max(double x_start,double x_end);
	double Min(double x_start,double x_end);
	double Peak();
	double ValleyLocal(double x_start,double x_end);
	double Valley();
	double PeakLocal(double x_start,double x_end);
	
	void Normalize(double new_max);
	void NormalizeByLocal(double x_start,double x_end,double new_max);
	void xShift(double dx);
	void Unwrap();
	void RemoveLinearTerm(double x_start,double x_end);

	int    FW(double& x1,double & x2,double threshold);
	double FW(double threshold);
	double FWStart(double threshold);
	double FWEnd(double threshold);
	double FWCenter(double threshold);
	double FWLocal(double x_start,double x_end,double threshold);
	double FWStartLocal(double x_start,double x_end,double threshold);
	double FWEndLocal(double x_start,double x_end,double threshold);
	double FWCenterLocal(double x_start,double x_end,double threshold);
	double FWHM();
	double FWHMCenter();
	double FWHMLocal(double x_start,double x_end);
	double TransitionPoint(double val);
	double TransitionPointLocal(double x_start,double x_end,double val);
	double TransitionInterval(double val1,double val2);
	double TransitionIntervalLocal(double x_start,double x_end,double val1,double val2);

	double Correl();
	double R2();
};

#endif // #ifndef CXYLIST_H_INCLUDED