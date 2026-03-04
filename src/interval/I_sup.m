function output = I_sup(var)
   global INTERVAL_MODE;
   if INTERVAL_MODE
      output = sup(var);
   else
      output = var;
   end
end



