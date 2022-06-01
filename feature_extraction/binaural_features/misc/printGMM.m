function printGMM(GMM, prefix)

for g=1:length(GMM)
    f = fopen(sprintf('%s_%d.txt', prefix, g), 'w');
    n = GMM{1}(1).nin;
    fprintf(f, 'n = %d\n', n);
    nCl = length(GMM{g});
    fprintf(f, 'nCl = %d\n', nCl);
    for cl=1:nCl
        nGauss = GMM{g}(cl).ncentres;
        fprintf(f, 'nGauss = %d\n', nGauss);
        for n=1:nGauss
            fprintf(f, 'prior = %.16e\n', GMM{g}(cl).priors(n));
            fprintf(f, 'center = ');
            printVect(f, GMM{g}(cl).centres(n, :));
            fprintf(f, 'covar_type = %s\n', GMM{g}(cl).covar_type);
            fprintf(f, 'covar = ');
            switch GMM{g}(cl).covar_type
                case 'spherical'
                    fprintf(f, '%.16e\n', GMM{g}(cl).covar(n));
                case 'diag'
                    printVect(f, GMM{g}(cl).covars(n, :));
                case 'full'
                    temp = GMM{g}(cl).covars(:, :, n);
                    printVect(f, temp(:));
                otherwise
                    error('Structure Error!');
            end
        end
    end
    fclose(f);
end

end

function printVect(f, x)

x = x(:);
fprintf(f, '[ ');
for n=1:length(x)
    fprintf(f, '%.16e ', x(n));
end
fprintf(f, ']\n');

end
            