
function comp = complement_set( set , space )

% comp is the set of things that are in space but not in set


comp = setxor( set , space );

end