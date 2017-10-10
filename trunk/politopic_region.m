function X = politopic_region(x, y, z)

X = Polyhedron('lb', [x(1), y(1), z(1)], 'ub', [x(2), y(2), z(2)]);

end

