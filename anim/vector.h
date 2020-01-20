struct vector
{
  double x;
  double y;
  double z;
  vector(double x, double y, double z) : x(x), y(y), z(z) {}
  vector() {}
};


  vector operator+(const vector& a, const vector& b)
  {
    return vector(a.x+b.x, a.y+b.y, a.z+b.z);;
  }

  double operator*(const vector& a, const vector& b)
  {
    return a.x*b.x + a.y*b.y + a.z*b.z;
  }

  vector operator^(const vector& a, const vector& b)
  {
    static vector c;
    c.x=a.y*b.z - a.z*b.y;
    c.y=a.z*b.x - a.x*b.z;
    c.z=a.x*b.y - a.y*b.x;
    return c;
  }

  vector operator*(const double& a, const vector& b)
  {
    static vector c;
    c.x = a*b.x;
    c.y = a*b.y;
    c.z = a*b.z;
    return c;
  }

  vector operator*(const vector& a, const double& b)
  {
    static vector c;
    c.x=a.x*b;
    c.y=a.y*b;
    c.z=a.z*b;
    return c;
  }

  vector operator-(const vector& a, const vector& b)
  {
    static vector c;
    c.x = a.x-b.x;
    c.y = a.y-b.y;
    c.z = a.z-b.z;
    return c;
  }

  vector operator/(const vector& a, const double& b)
  {
    static vector c;
    c.x = a.x/b;
    c.y = a.y/b;
    c.z = a.z/b;
    return c;
  }

  vector operator+=(vector& a, const vector b)
    {
    a=a+b;
    return a;
    }

  vector operator-=(vector& a, const vector b)
    {
    a=a-b;
    return a;
    }

  double radius( vector BD1, vector BD2 ){
    return sqrt( (BD1 - BD2) * (BD1 - BD2) );//FIXME Does this file need to include math.h to get sqrt()?
  }

  double mag (vector v)
  {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
  }
