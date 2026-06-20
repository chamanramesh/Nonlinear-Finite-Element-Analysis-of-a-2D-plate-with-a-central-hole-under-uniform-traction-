function[A_e] = Assembly(element,A_e)

for j = 1:numel(element)
    node_id = element(j);
    u = 2 * node_id - 1;
    v = 2 * node_id;
    A_e(j * 2 - 1, u) = 1;
    A_e(j * 2, v) = 1;
end
end
