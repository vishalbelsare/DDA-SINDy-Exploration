function [value, isterminal, direction] = events_diverge(t,y);

value = max(abs(y))-1e5;
isterminal = 1;
direction = [];

