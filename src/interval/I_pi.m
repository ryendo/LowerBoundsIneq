function output = I_pi
   global INTERVAL_MODE;
   if INTERVAL_MODE
      % Use the string form so the interval rigorously encloses the true
      % constant pi. NOTE: intval(pi) would enclose only the (rounded)
      % double-precision value of pi, which does NOT contain the exact pi.
      output = intval('pi');
   else
      output = pi;
   end
end



