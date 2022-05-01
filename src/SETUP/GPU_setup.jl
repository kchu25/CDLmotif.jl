# kernel thread-block setup
const threads_1d = 512;
const threads_2d = 32;
const threads_3d = 10;
const ker_1d = threads_1d;
const ker_2d = (threads_2d, threads_2d);
const ker_3d = (threads_3d, threads_3d, threads_3d);

b_size_1d(X) = ceil.(Int, size(X) ./ threads_1d)
b_size_2d(X) = ceil.(Int, size(X) ./ threads_2d)
b_size_3d(X) = ceil.(Int, size(X) ./ threads_3d)