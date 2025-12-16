function [output] = derivativeDet(M,Mi,option)


if option == "inv"

    output = trace(M*Mi);

else

    output = trace(M\Mi);

end

end