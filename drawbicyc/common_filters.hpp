#pragma once

using namespace blitz;

double iscloudy_rc(double x)
{
  return x > 1.e-4 ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(iscloudy_rc)

double is_th_prtrb(double x)
{
  return x > 300.1 ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(is_th_prtrb)

double iscloudy(double x)
{
  return x > 20. ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(iscloudy)

double isdowndraught(double x)
{
  return  x < -0.2 ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(isdowndraught)

double isupdraught(double x)
{
  return  x > 0.2 ? 1. : 0.; 
}
BZ_DECLARE_FUNCTION(isupdraught)
