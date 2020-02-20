function [Kripple] = filter_images(Knormalized)
    Kripple = zeros(size(Knormalized));
    p = Progress(size(Knormalized,3));
    for i = 1:size(Knormalized,3)
        p.d(i)
        A = Knormalized(:,:,i);
        A = A - imfilter(double(A),fspecial('average',100),'replicate');
        A = imfilter(double(A),fspecial('gaussian',100,25),'replicate');
        A = A - mean(A(:));
        Kripple(:,:,i) = A;
    end
end