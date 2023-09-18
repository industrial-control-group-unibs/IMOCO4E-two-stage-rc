function friction=LUTfriction(position,velocity,LUT_position,LUT_Coloumb_friction,LUT_viscous_friction)

if (position<=LUT_position(1))
    Coloumb_friction=LUT_Coloumb_friction(1);
    viscous_friction=LUT_viscous_friction(1);
elseif (position>=LUT_position(end))
    Coloumb_friction=LUT_Coloumb_friction(end);
    viscous_friction=LUT_viscous_friction(end);
else
    for idx=2:length(LUT_position)
        if (position<=LUT_position(idx))
            x=(position-LUT_position(idx-1))/(LUT_position(idx)-LUT_position(idx-1));
            Coloumb_friction=LUT_Coloumb_friction(idx-1)+(LUT_Coloumb_friction(idx)-LUT_Coloumb_friction(idx-1))*x;
            viscous_friction=LUT_viscous_friction(idx-1)+(LUT_viscous_friction(idx)-LUT_viscous_friction(idx-1))*x;
            break;
        end
    end
end

direction=tanh(velocity*100.0);  % approximation of sign

friction=Coloumb_friction*direction+viscous_friction*velocity;
