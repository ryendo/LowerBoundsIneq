function output = I_intval(var)
   global INTERVAL_MODE;
   % Debug: check INTERVAL_MODE value
   if ~exist('INTERVAL_MODE', 'var') || isempty(INTERVAL_MODE)
       fprintf('[I_intval DEBUG] INTERVAL_MODE is not set or empty\n');
       INTERVAL_MODE = 0;
   end
   if INTERVAL_MODE
      output = intval(var);
   else
      if ischar(var)
        output = str2double(var);
      else
        output = var;
      end
   end
end



