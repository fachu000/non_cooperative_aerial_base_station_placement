function str = print_frequency(f)

if (0<=f)&&(f<1e3)
	str = sprintf('%g Hz',f);
elseif (1e3<=f)&&(f<1e6)
	str = sprintf('%g KHz',f/1e3);
elseif (1e6<=f)&&(f<1e9)
	str = sprintf('%g MHz',f/1e6);
elseif (1e9<=f)
	str = sprintf('%g GHz',f/1e9);
end

end