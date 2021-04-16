function [prob_inclusion, alpha, beta, gamma, beta0, sigmasq, output] = ASI_sampler_special(data_fixed, data, target, g, gprior, mode, RB_adap, RB, hparam, tau, nu, burnin, numbofits, thin, numbofreps, heat, adap_type, updateg, fig1)

    warning off
    
    if ( verLessThan('matlab', '9.7') == 0 )
        d = uiprogressdlg(fig1,'Title','Fitting Model','Message', 'Progress');
    else
        uialert(fig1, 'Running Code', 'Fitting Run','Icon', 'info');
    end
    
    % disp(['hparam = ' num2str(hparam)]);
    % disp(['g = ' num2str(g)]);
    
    rhostar = zeros(1, length(heat) - 1);
    for i = 1:(length(heat)-1)
        rhostar(i) = log(log(heat(i)/heat(i+1)));
    end
    order = zeros(numbofreps, length(heat));
    for rep = 1:numbofreps
        order(rep, :) = 1:length(heat);
    end
    
    numb_subjects = size(data, 1);
    n = zeros(1, numb_subjects);
    p = zeros(1, numb_subjects);
    for subj = 1:numb_subjects
        n(subj) = length(target{subj, 1});
        p(subj) = size(data{subj, 1}, 2);
    end
    
    numbofchains = length(heat);
    
    data_fixed2 = zeros(sum(n), size(data_fixed{1, 1}, 2));
    target2 = zeros(sum(n), 1);
    count = 0;
    for subj = 1:numb_subjects    
        data_fixed2((count + 1):(count + n(subj)), :) = data_fixed{subj, 1};
        target2((count + 1):(count + n(subj)), :) = target{subj, 1};
        
        count = count + n(subj);
    end
    
    % F = eye(sum(n)) - data_fixed2 * inv(data_fixed2' * data_fixed2) * data_fixed2';
    % target2star = F * target2;
    
    alpha_fixed = zeros(size(data_fixed{1, 1}, 2), numbofreps, numbofchains);
    resid = cell(numb_subjects, 1);
    
    % mat1 = inv(data_fixed2' * data_fixed2);
    % for subj = 1:numb_subjects
    %     data{subj, 1} = data{sub, j} - data_fixed{subj, 1} * mat1 * (data_fixed2' * data);
    % target = target - data_fixed2 * mat1 * (data_fixed2' * target2);
    
    
    
    for rep = 1:numbofreps
        for chain = 1:numbofchains
            alpha_fixed(:, rep, chain) = inv(data_fixed2' * data_fixed2) * (data_fixed2' * target2);
        end
    end
    
    for subj = 1:numb_subjects
        resid{subj, 1} = zeros(n(subj), numbofreps, numbofchains);
        for rep = 1:numbofreps
            for chain = 1:numbofchains
                resid{subj, 1}(:, rep, chain) = target{subj, 1} - data_fixed{subj, 1} * alpha_fixed(:, rep, chain);
            end
        end
    end
    
    
    % beta_sat = regress(target, [ones(n, 1) data]);
    % RSS = sum((target - [ones(n, 1) data] * beta_sat).^2);
    % sigmasqhat1 = RSS / size(data, 1);
    % lambda = sigmasqhat1;
    % g = 10 / lambda;
    
    % sigmasqalpha = 4 / 2;
    % sigmasqbeta = 4 / 2 * lambda;
    sigmasqalpha = 10 * ones(numbofreps, numbofchains);
    sigmasqbeta = 10 * ones(numbofreps, numbofchains);
    
    g = g * ones(numbofreps, numbofchains);
    g_fixed = g;
    sigmasq_int = ones(numbofreps, numbofchains);
    
    if ( length(hparam) == 1 )
        fixed = 1;
        w = hparam * ones(numbofreps, numbofchains);
        wstar = hparam;
    else
        fixed = 0;
        wa = hparam(1);
        wb = hparam(2);
        wstar = wa / (wa + wb);
    end
    
    fixed
    
    logita = zeros(numb_subjects, 1);
    logitb = zeros(numb_subjects, 1);
    for subj = 1:numb_subjects
        logita(subj, 1) = 0.1 / p(subj);
        logitb(subj, 1) = 1 - 0.1 / p(subj);
    end
    
    % XTy_all = cell(numb_subjects, 1);
    % for subj = 1:numb_subjects
    %     XTy_all{subj, 1} = data{subj, 1}' * target{subj, 1};
    % end
    
    loglike = zeros(numb_subjects, numbofreps, numbofchains);
    
    C = cell(numb_subjects, numbofreps, numbofchains);
    gamma = cell(numb_subjects, 1);
    for subj = 1:numb_subjects
        gamma{subj, 1} = zeros(numbofreps, p(subj), numbofchains);
        for chain = 1:numbofchains
            for rep = 1:numbofreps
                check = 0;
                while ( check == 0 )
                    gamma{subj, 1}(rep, :, chain) = rand(1, p(subj)) < wstar;
                    
                    datastar = [ones(n(subj), 1) data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)];
                    if ( gprior == 1 )
                        n0star = zeros(sum(gamma{subj, 1}(rep, :, chain)));
                        n0star((1 + size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma{subj, 1}(rep, :, chain))), (1 + size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma{subj, 1}(rep, :, chain)))) = 1 / g(rep, chain)  * datastar(:, 2:end)' * datastar(:, 2:end);
                    else
                        n0star = 1 / g(rep, chain) * eye(sum(gamma{subj, 1}(rep, :, chain)) + 1);
                        n0star(1, 1) = 1 / sigmasq_int(rep, chain);
                    end
                    C{subj, rep, chain} = inv(datastar' * datastar + n0star);
                    
                    
                    XTy = [sum(resid{subj, 1}(:, rep, chain)); data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)' * resid{subj, 1}(:, rep, chain)];
                    loglike(subj, rep, chain) = - 0.5 * sum(gamma{subj, 1}(rep, :, chain)) * log(g(rep, chain));
                    loglike(subj, rep, chain) = loglike(subj, rep, chain) - 0.5 * log(det(datastar' * datastar + n0star));
                    loglike(subj, rep, chain) = loglike(subj, rep, chain) - (0.5 * n(subj) + sigmasqalpha(rep, chain)) * log(sigmasqbeta(rep, chain) + 0.5 * (resid{subj, 1}(:, rep, chain)' * resid{subj, 1}(:, rep, chain)) - 0.5 * XTy' * C{subj, rep, chain} * XTy);
                    
                    if ( isreal(loglike(subj, rep, chain)) == 1 )
                        check = 1;
                    end
                end
                if ( size(C{subj, rep, chain}, 2) ~=  1 + sum(gamma{subj, 1}(rep, :, chain)))
                    size(C{subj, rep, chain})
                    size(datastar)
                    size(n0star)
                    sum(gamma{subj, 1}(rep, :, chain))
                    stop;
                end
    
            end
        end
    end
    
    logitzeta0 = zeros(numb_subjects, numbofchains);
    for subj = 1:numb_subjects
        zetahat = nu / (wstar * p(subj));
        if ( zetahat > 1 )
            zetahat = 0.1;
        end
        logitzeta0(subj, :) = (log(zetahat - logita(subj)) - log(logitb(subj) - zetahat)) * ones(1, numbofchains);
    end
    totalaccept = zeros(1, numbofchains);
    totalcount = zeros(1, numbofchains);
    totalchange = zeros(1, numbofchains);
    totalchange2 = zeros(1, numbofchains);
    totalaccept2 = zeros(1, numbofchains - 1);
    totalcount2 = zeros(1, numbofchains - 1);
    gaccept = zeros(1, numbofchains);
    gcount = zeros(1, numbofchains);
    g_fixedaccept = zeros(1, numbofchains);
    g_fixedcount = zeros(1, numbofchains);
    sigmasqalphaaccept = zeros(1, numbofchains);
    sigmasqalphacount = zeros(1, numbofchains);
    
    loggsd = zeros(1, numbofchains);
    logg_fixedsd = zeros(1, numbofchains);
    logsigmasqalphasd = log(0.01) * ones(1, numbofchains);
    
    sumgamma = cell(numb_subjects, 1);
    sumgamma2 = cell(numb_subjects, 1);
    sumgamma3 = cell(numb_subjects, 1);
    count3 = cell(numb_subjects, 1);
    for subj = 1:numb_subjects
        sumgamma{subj, 1} = zeros(1, p(subj));
        sumgamma2{subj, 1} = zeros(1, p(subj));
        sumgamma3{subj, 1} = zeros(numbofchains, p(subj));
        count3{subj, 1} = zeros(numbofchains, 1);
    end
    %holdloglike = zeros(numb_subjects, numbofits * numbofreps);
    %holdaccept = zeros(numb_subjects, numbofits);
    holdg = zeros(1, numbofits * numbofreps);
    holdg_fixed = zeros(1, numbofits * numbofreps);
    holdalpha = zeros(numb_subjects, numbofreps * numbofits);
    if ( mode == 2 )
        holdgamma = cell(numb_subjects, numbofreps * numbofits);
        holdbeta = cell(numb_subjects, numbofreps * numbofits);
    else
        holdalpha = [];
        holdbeta = [];
    end
    holdbeta0 = zeros(numbofreps * numbofits, size(data_fixed{1, 1}, 2));
    holdsigmasqalpha = zeros(1, numbofreps * numbofits);
    holdsigmasqbeta = zeros(1, numbofreps * numbofits);
    holdsigmasq = zeros(numb_subjects, numbofreps * numbofits);
    holdsigmasq_int = zeros(1, numbofreps * numbofits);
    
    for it = 1:(burnin+numbofits*thin)
        
        if ( mod(it, 10) == 0 )
            top_chain = find(order(rep, :) == 1);
    
            disp(['it = ' num2str(it)]);
            disp(['numbofreps = ' num2str(numbofreps)]);
            disp(['logitzeta0 = ' num2str(logitzeta0(chain1))]);
            disp(['accept = ' num2str(totalaccept./totalcount)]);
            disp(['accept* = ' num2str(totalaccept./totalchange2)]);
            disp(['g accept = ' num2str(gaccept./gcount)]);
            disp(['sigmasqalpha accept = ' num2str(sigmasqalphaaccept./sigmasqalphacount)]);
            disp(['sigmasqalpha = ' num2str(sigmasqalpha(:, top_chain)')]);
            disp(['sigmasqbeta = ' num2str(sigmasqbeta(:, top_chain)')]);
            disp(['g = ' num2str(g(:, top_chain)')]);
            disp(['tau = ' num2str(tau)]);
            disp(['change = ' num2str(totalchange./totalcount)]);
            disp(num2str(heat));
            disp(num2str(totalaccept2./totalcount2));
            disp(' ');
        end
        
        for chain = 1:numbofchains
            for rep = 1:numbofreps
                
                chain1 = order(rep, chain);
                
                
                sigmasq_all = zeros(1, numb_subjects);
                for subj = 1:numb_subjects
                    XTy = [sum(resid{subj, 1}(:, rep, chain)); data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)' * resid{subj, 1}(:, rep, chain)];
                    
                    alphastar = sigmasqalpha(rep, chain) + 0.5 * n(subj);
                    betastar = sigmasqbeta(rep, chain) + 0.5 * (resid{subj, 1}(:, rep, chain)' * resid{subj, 1}(:, rep, chain) - XTy' * C{subj, rep, chain} * XTy);
                    sigmasq_all(subj) = 1 / gamrnd(alphastar, 1 / betastar);
                    
                    if ( (isnan(sigmasq_all(subj)) == 1) || (isinf(sigmasq_all(subj)) == 1) )
                        alphastar
                        betastar
                        stop;
                    end
                end
                
                a1 = sigmasqalpha(rep, chain) * numb_subjects;
                b1 = sum(1./sigmasq_all);
                
                sigmasqbeta(rep, chain) = gamrnd(a1, 1 / b1);
                
                newsigmasqalpha = sigmasqalpha(rep, chain) * exp(exp(logsigmasqalphasd(chain1)) * randn);
                
                loglikestar = numb_subjects * (sigmasqalpha(rep, chain) * log(sigmasqbeta(rep, chain)) - gammaln(sigmasqalpha(rep, chain)));
                loglikestar = loglikestar + (sigmasqalpha(rep, chain) - 1) * sum(log(1 ./ sigmasq_all));
                
                newloglikestar = numb_subjects * (newsigmasqalpha * log(sigmasqbeta(rep, chain)) - gammaln(newsigmasqalpha));
                newloglikestar = newloglikestar + (newsigmasqalpha - 1) * sum(log(1 ./ sigmasq_all));
                
                logaccept = newloglikestar - loglikestar;
                
                accept = 1;
                if ( isreal(logaccept) == 0 )
                    accept = 0;
                elseif ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
                    accept = 0;
                elseif ( logaccept < 0 )
                    accept = exp(logaccept);
                end
                
                if ( rand < accept )
                    sigmasqalpha(rep, chain) = newsigmasqalpha;
                end
                
                sigmasqalphaaccept(chain1) = sigmasqalphaaccept(chain1) + accept;
                sigmasqalphacount(chain1) = sigmasqalphacount(chain1) + 1;
                
                logsigmasqalphasd(chain1) = logsigmasqalphasd(chain1) + 1 / it^0.55 * (accept - 0.3);
                
                
                
                
                
                for subj = 1:numb_subjects
                    if ( size(C{subj, rep, chain}, 2) ~=  1 + sum(gamma{subj, 1}(rep, :, chain)))
                        stop;
                    end
                end
                
                
                part1 = zeros(size(data_fixed2, 2), size(data_fixed2, 2));
                part2 = zeros(size(data_fixed2, 2), 1);
                for subj = 1:numb_subjects
                    datastar = [ones(n(subj), 1) data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)];
                    
                    BZ = data_fixed{subj, 1} - datastar *  C{subj, rep, chain} * (datastar' * data_fixed{subj, 1});
                    By = target{subj, 1} - datastar *  C{subj, rep, chain} * (datastar' * target{subj, 1});
                    part1 = part1 + data_fixed{subj, 1}' * data_fixed{subj, 1} / sigmasq_all(subj);
                    part1 = part1 - (data_fixed{subj, 1}' * datastar) * C{subj, rep, chain} * (datastar' * data_fixed{subj, 1}) / sigmasq_all(subj);
                    part2 = part2 + data_fixed{subj, 1}' * target{subj, 1} / sigmasq_all(subj);
                    part2 = part2 - (data_fixed{subj, 1}' * datastar) * C{subj, rep, chain} * (datastar' * target{subj, 1}) / sigmasq_all(subj);
                end
                %            varstar = inv(part1 + 1 / g_fixed(rep, chain) * eye(size(data_fixed{1, 1}, 2)));
                varstar = inv(part1);
                mustar = varstar * part2;
                alpha_fixed(:, rep, chain) = mustar + chol(varstar)' * randn(size(part1, 1), 1);
                            
                for subj = 1:numb_subjects
                    resid{subj, 1}(:, rep, chain) = target{subj, 1} - data_fixed{subj, 1} * alpha_fixed(:, rep, chain);
                end
                
                nu1 = 1 / gamrnd(1, 1 / (1 + 1 / sigmasq_int(rep, chain)));
                
                alphastar = 0.5 + 0.5 * numb_subjects;
                betastar = 1 / nu1;
                for subj = 1:numb_subjects
                    XTy = [sum(resid{subj, 1}(:, rep, chain)); data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)' * resid{subj, 1}(:, rep, chain)];
                    
                    mustar = C{subj, rep, chain} * XTy;
                    varstar = sigmasq_all(subj) * C{subj, rep, chain};
                    
                    beta_local = mustar + chol(varstar)' * randn(size(C{subj, rep, chain}, 1), 1);                
    
                    betastar = betastar + 0.5 * beta_local(1)^2 / sigmasq_all(subj);
                end
                
                sigmasq_int(rep, chain) = 1 / gamrnd(alphastar, 1 / betastar);
    
                
                
                %             newg_fixed = g_fixed(rep, chain) * exp(exp(loggsd_fixed(chain1)) * randn);
                %
                %             loglikestar = - 0.5 * (size(data_fixed{1, 1}, 2) - 1) * log(g_fixed(rep, chain)) - 0.5 * sum(alpha_fixed(2:end, rep, chain).^2) / g_fixed(rep, chain);
                %             newloglikestar = - 0.5 * (size(data_fixed{1, 1}, 2) - 1) * log(newg_fixed) - 0.5 * sum(alpha_fixed(2:end, rep, chain).^2) / newg_fixed;
                %
                %             logaccept =  newloglikestar - loglikestar;
                %
                %             logaccept = logaccept + 0.5 * log(newg_fixed) - log(1 + newg_fixed);
                %             logaccept = logaccept - 0.5 * log(g_fixed(rep, chain)) + log(1 + g_fixed(rep, chain));
                %
                %             accept = 1;
                %             if ( isreal(logaccept) == 0 )
                %                 accept = 0;
                %             elseif ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
                %                 accept = 0;
                %             elseif ( logaccept < 0 )
                %                 accept = exp(logaccept);
                %             end
                %
                %             if ( rand < accept )
                %                 g_fixed(rep, chain) = newg_fixed;
                %             end
                %
                %             g_fixedaccept(chain1) = g_fixedaccept(chain1) + accept;
                %             g_fixedcount(chain1) = g_fixedcount(chain1) + 1;
                %
                %             logg_fixedsd(chain1) = logg_fixedsd(chain1) + 1 / it^0.55 * (accept - 0.234);
                
                
                
                for subj = 1:numb_subjects
                    
                    datastar = [ones(n(subj), 1) data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)];
                    if ( gprior == 1 )
                        n0star = zeros(sum(gamma(rep, :, chain)));
                        n0star((1+size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma(rep, :, chain))), (1+size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma(rep, :, chain)))) = 1 / g(rep, chain) * datastar(:, 2:end)' * datastar(:, 2:end);
                    else
                        n0star = 1 / g(rep, chain) * eye(sum(gamma{subj, 1}(rep, :, chain)) + 1);
                        n0star(1, 1) = 1 / sigmasq_int(rep, chain);
                    end
                    C{subj, rep, chain} = inv(datastar' * datastar + n0star);
    
                    XTy = [sum(resid{subj, 1}(:, rep, chain)); data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)' * resid{subj, 1}(:, rep, chain)];
                    loglike(subj, rep, chain) = - 0.5 * sum(gamma{subj, 1}(rep, :, chain)) * log(g(rep, chain));
                    loglike(subj, rep, chain) = loglike(subj, rep, chain) - 0.5 * log(det(datastar' * datastar + n0star));
                    loglike(subj, rep, chain) = loglike(subj, rep, chain) - (0.5 * n(subj) + sigmasqalpha(rep, chain)) * log(sigmasqbeta(rep, chain) + 0.5 * (resid{subj, 1}(:, rep, chain)' * resid{subj, 1}(:, rep, chain)) - 0.5 * XTy' * C{subj, rep, chain} * XTy);
                    
                                  
                    if ( logitzeta0(subj, chain1) < 0 )
                        zeta = (logita(subj) + logitb(subj) * exp(logitzeta0(subj, chain1))) ./ (1 + exp(logitzeta0(subj, chain1)));
                    else
                        zeta = (logita(subj) * exp(- logitzeta0(subj, chain1)) + logitb(subj)) ./ (exp(- logitzeta0(subj, chain1)) + 1);
                    end
                    probstar = sumgamma3{subj, 1}(chain1, :) / count3{subj, 1}(chain1);
                    probstar = 0.0001 + (1 - 0.0002) * probstar;
                    
                    p1 = probstar;
                    p2 = 1 - probstar;
                    
                    zetaAS = p1 .* zeta .* min(1./p1, 1./p2);
                    zetaDS = p2 .* zeta .* min(1./p1, 1./p2);
                    
                    gammaA = zeros(1, p(subj));
                    gammaD = zeros(1, p(subj));
                    newgamma = gamma{subj, 1}(rep, :, chain);
                    logprob = 0;
                    
                    gammaA(gamma{subj, 1}(rep, :, chain)==0) = rand(1, sum(gamma{subj, 1}(rep, :, chain)==0)) < zetaAS(gamma{subj, 1}(rep, :, chain)==0);
                    newgamma(gammaA==1) = 1;
                    gammaD(gamma{subj, 1}(rep, :, chain)==1) = rand(1, sum(gamma{subj, 1}(rep, :, chain)==1)) < zetaDS(gamma{subj, 1}(rep, :, chain)==1);
                    newgamma(gammaD==1) = 0;
                    temp = (gamma{subj, 1}(rep, :, chain)==0) .* (newgamma==1);
                    logprob = logprob - sum(log(zetaAS(temp==1)));
                    temp = (gamma{subj, 1}(rep, :, chain)==0) .* (newgamma==0);
                    logprob = logprob - sum(log(1 - zetaAS(temp==1)));
                    temp = (gamma{subj, 1}(rep, :, chain)==1) .* (newgamma==0);
                    logprob = logprob - sum(log(zetaDS(temp==1)));
                    temp = (gamma{subj, 1}(rep, :, chain)==1) .* (newgamma==1);
                    logprob = logprob - sum(log(1 - zetaDS(temp==1)));
                    temp = (newgamma==0) .* (gamma{subj, 1}(rep, :, chain)==1);
                    logprob = logprob + sum(log(zetaAS(temp==1)));
                    temp = (newgamma==0) .* (gamma{subj, 1}(rep, :, chain)==0);
                    logprob = logprob + sum(log(1 - zetaAS(temp==1)));
                    temp = (newgamma==1) .* (gamma{subj, 1}(rep, :, chain)==0);
                    logprob = logprob + sum(log(zetaDS(temp==1)));
                    temp = (newgamma==1) .* (gamma{subj, 1}(rep, :, chain)==1);
                    logprob = logprob + sum(log(1 - zetaDS(temp==1)));
                    change = sum(gammaA) + sum(gammaD);
                    
                    if ( change == 0 )
                        newloglike = loglike(subj, rep, chain);
                    else
                        datastar = [ones(n(subj), 1) data{subj, 1}(:, newgamma==1)];
                        if ( gprior == 1 )
                            n0star = zeros(sum(newgamma));
                            n0star((1 + size(data_fixed, 2)):(sum(newgamma) + size(data_fixed, 2)), (1 + size(data_fixed, 2)):(sum(newgamma) + size(data_fixed, 2))) = 1 / g(rep, chain) * datastar(:, 2:end)' * datastar(:, 2:end);
                        else
                            n0star = 1 / g(rep, chain) * eye(sum(newgamma) + 1);
                            n0star(1, 1) = 1 / sigmasq_int(rep, chain);
                        end
                        newC = inv(datastar' * datastar + n0star);
                        
                        XTy = [sum(resid{subj, 1}(:, rep, chain)); data{subj, 1}(:, newgamma==1)' * resid{subj, 1}(:, rep, chain)];
                        newloglike = - 0.5 * sum(newgamma) * log(g(rep, chain));
                        newloglike = newloglike - 0.5 * log(det(datastar' * datastar + n0star));
                        newloglike = newloglike - (0.5 * n(subj) + sigmasqalpha(rep, chain)) * log(sigmasqbeta(rep, chain) + 0.5 * (resid{subj, 1}(:, rep, chain)' * resid{subj, 1}(:, rep, chain)) - 0.5 * XTy' * newC * XTy);
                    end
                    
                    logaccept = heat(chain1) * (newloglike - loglike(subj, rep, chain));
                    if ( fixed == 1 )
                        logaccept = logaccept + sum(newgamma) * log(w(rep, chain)) + (p - sum(newgamma)) * log(1 - w(rep, chain));
                        logaccept = logaccept - sum(gamma{subj, 1}(rep, :, chain)) * log(w(rep, chain)) - (p - sum(gamma{subj, 1}(rep, :, chain))) * log(1 - w(rep, chain));
                    else
                        logaccept = logaccept + gammaln(wa + sum(newgamma)) + gammaln(wb + p(subj) - sum(newgamma));
                        logaccept = logaccept - gammaln(wa + sum(gamma{subj, 1}(rep, :, chain))) - gammaln(wb + p(subj) - sum(gamma{subj, 1}(rep, :, chain)));
                    end
                    
                    logaccept = logaccept + logprob;
                    
                    accept = 1;
                    if ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) || (isreal(logaccept) == 0) )
                        accept = 0;
                    elseif ( logaccept < 0 )
                        accept = exp(logaccept);
                    end
                    
                    if ( change ~= 0 )
                        if ( rand < accept )
                            gamma{subj, 1}(rep, :, chain) = newgamma;
                            loglike(subj, rep, chain) = newloglike;
                            C{subj, rep, chain} = newC;
                        end
                    end
                    
                    totalaccept(chain1) = totalaccept(chain1) + (change~=0)*accept;
                    totalchange(chain1) = totalchange(chain1) + change;
                    totalchange2(chain1) = totalchange2(chain1) + (change~=0);
                    totalcount(chain1) = totalcount(chain1) + 1;
                    
                    %                             for rep1 = 1:numbofreps
                    %                 for chain1 = 1:numbofchains
                    %                     for subj1 = 1:numb_subjects
                    %                         if ( size(C{subj1, rep1, chain1}, 2) ~=  1 + sum(gamma{subj1, 1}(rep1, :, chain1)))
                    %                             rep1
                    %                             chain1
                    %                             subj1
                    %                             stop;
                    %                         end
                    %                     end
                    %                 end
                    %             end
                    %
                    %
                    if ( (adap_type == 1) || ((adap_type==0) && (it < burnin)) )
                        
                        logitzeta0(subj, chain1) = logitzeta0(subj, chain1) + 2 / it^0.55 * (change~=0) * (accept - tau);
                        
                        if ( isreal(logitzeta0(subj, chain1)) == 0 )
                            accept
                            logitzeta0(chain1)
                            stop;
                        end
                        
                        if ( (it > burnin)  )
                            sumgamma3{subj, 1}(chain1, :) = sumgamma3{subj, 1}(chain1, :) + gamma{subj, 1}(rep, :, chain);
                            count3{subj, 1}(chain1) = count3{subj, 1}(chain1) + 1;
                        else
                            datastar = [ones(n(subj), 1) data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)];
                            if ( gprior == 1 )
                                n0star = zeros(sum(gamma{subj, 1}(rep, :, chain)));
                                n0star((1+size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma(rep, :, chain))), (1+size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma(rep, :, chain)))) = 1 / g(rep, chain) * datastar(:, 2:end)' * datastar(:, 2:end);
                            else
                                n0star = 1 / g(rep, chain) * eye(sum(gamma{subj, 1}(rep, :, chain)) + 1);
                                n0star(1, 1) = 1 / sigmasq_int(rep, chain);
                            end
                            
                            XTy = [sum(resid{subj, 1}(:, rep, chain)); data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)' * resid{subj, 1}(:, rep, chain)];
                            inner1 = XTy' * C{subj, rep, chain} * XTy;
                            
                            oldloglike2 = 0;
                            oldloglike2 = oldloglike2 - (0.5 * n(subj) + sigmasqalpha(rep, chain)) * log(sigmasqbeta(rep, chain) + 0.5 * (resid{subj, 1}(:, rep, chain)' * resid{subj, 1}(:, rep, chain)) - 0.5 * inner1);
                            
                            if ( isreal(oldloglike2) == 0 )
                                resid' * resid
                                XTy'
                                g
                                C{rep, chain}
                                inner1
                                stop;
                            end
                            
                            zstar = find(gamma{subj, 1}(rep, :, chain) == 1 );
                            
                            XTx_all = datastar' * data{subj, 1};
                            XTy_all = data{subj, 1}' * resid{subj, 1}(:, rep, chain);
                            
                            D = zeros(1, p(subj));
                            newloglike2 = zeros(1, p(subj));
                            
                            Gstar = XTy' * C{subj, rep, chain} * XTx_all;
                            Gstar2 = 2 * resid{subj, 1}(:, rep, chain)' * data{subj, 1};
                            %                     Astar = XTy_all .* Gstar';
                            Cstar = XTy_all .* XTy_all;
                            
                            Bstar = zeros(1, p(subj));
                            for j = 1:p(subj)
                                XTx = XTx_all(:, j);
                                
                                if ( gamma{subj, 1}(rep, j, chain) == 0 )
                                    D(j) = sum(data{subj, 1}(:, j).^2) + 1 / g(rep, chain) - XTx' * C{subj, rep, chain} * XTx;
                                end
                            end
                            
                            x1 = gamma{subj, 1}(rep, :, chain)==0;
                            %                    inner2 = (Gstar(x1==1) - 2 * XTy_all(x1==1)') .* Gstar(x1==1) + Cstar(x1==1)';
                            inner2 = (Gstar(x1==1) - Gstar2(x1==1)) .* Gstar(x1==1) + Cstar(x1==1)';
                            %                    inner2 = inner2 - 2 * Astar(x1==1)';
                            %                   inner2 = inner2 - 2 * XTy_all(x1 == 1)' .* Gstar(x1==1);
                            %              inner2 = inner2 + Cstar(x1==1)';
                            Bstar(x1==1) = inner1 + inner2 ./ D(x1==1);
                            newloglike2(x1==1) = - (0.5 * n(subj) + sigmasqalpha(rep, chain)) * log(sigmasqbeta(rep, chain) + 0.5 * (resid{subj, 1}(:, rep, chain)' * resid{subj, 1}(:, rep, chain)) - 0.5 * Bstar(x1==1));
                            
                            if ( isreal(newloglike2) == 0 )
                                D(x1==1)
                                target' * target - Bstar(x1==1)
                                Bstar(x1==1)
                                newloglike2(x1==1)
                                stop;
                            end
                            
                            
                            for j = 1:p
                                XTx = XTx_all(:, j);
                                
                                if ( gamma{subj, 1}(rep, j, chain) == 1 )
                                    pos = find(zstar==j) + 1;
                                    Fstar = C{subj, rep, chain}([1:(pos-1) (pos+1):end], [1:(pos-1) (pos+1):end]);
                                    d3 = C{subj, rep, chain}(pos, pos);
                                    u3 = - C{subj, rep, chain}([1:(pos-1) (pos+1):end], pos);
                                    %                            u2 = u3 / d3;
                                    %                            inv2 = Fstar - d3 * u2 * u2';
                                    inv2 = Fstar - (u3 * u3') / d3;
                                    
                                    D(j) = sum(data{subj, 1}(:, j).^2) + 1 / g(rep, chain) - XTx([1:(pos-1) (pos+1):end], 1)' * inv2 * XTx([1:(pos-1) (pos+1):end], 1);
                                    
                                    newXTy = [XTy(1:(pos-1), 1); XTy((pos+1):end, 1)];
                                    newloglike2(j) = - (0.5 * n(subj) + sigmasqalpha(rep, chain)) * log(sigmasqbeta(rep, chain) + 0.5 * (resid{subj, 1}(:, rep, chain)' * resid{subj, 1}(:, rep, chain)) - 0.5 * newXTy' * inv2 * newXTy);
                                end
                                
                                if ( isreal(newloglike2(j)) == 0 )
                                    target{subj, 1}' * target{subj, 1} - 0.5 * newXTy' * inv2 * newXTy
                                    j
                                    newloglike2(j)
                                    D(j)
                                end
                            end
                            
                            if ( isreal(oldloglike2) == 0 )
                                subj
                                stop;
                            end
                            
                            if ( isreal(newloglike2) == 0 )
                                subj
                                inv2
                                stop;
                            end
                            
                            newloglike2 = newloglike2 - (2 * (gamma{subj, 1}(rep, :, chain)==0) - 1) .* (0.5 * log(g(rep, chain)) + 0.5 * log(abs(D)));
                            %                    newloglike2 = newloglike2 - (2 * (gamma(rep, :, chain)==0) - 1) .* (0.5 * log(g(rep, chain)) + 0.5 * log(D));
                            
                            if ( isreal(newloglike2) == 0 )
                                g
                                D
                                stop;
                            end
                            
                            diff = heat(chain1) * (2 * (gamma{subj, 1}(rep, :, chain)==0) - 1) .* (newloglike2 - oldloglike2);
                            
                            if ( isreal(diff) == 0 )
                                stop;
                            end
                            
                            BF = zeros(1, p(subj));
                            if ( fixed == 1 )
                                x1 = (diff < 0);
                                if ( sum(x1) > 0 )
                                    BF(x1==1) = w(rep, chain) * exp(diff(x1==1)) ./ (w(rep, chain) * exp(diff(x1==1)) + 1 - w(rep, chain));
                                end
                                
                                x1 = (diff >= 0);
                                if ( sum(x1) > 0 )
                                    BF(x1==1) = w(rep, chain) ./ (w(rep, chain) + (1 - w(rep, chain)) * exp(- diff(x1==1)));
                                end
                            else
                                wstar = (wa + sum(gamma{subj, 1}(rep, :, chain))) / (wa + wb + p(subj));
                                x1 = (diff < 0);
                                if ( sum(x1) > 0 )
                                    BF(x1==1) = wstar * exp(diff(x1==1)) ./ (wstar * exp(diff(x1==1)) + 1 - wstar);
                                end
                                
                                x1 = (diff >= 0);
                                if ( sum(x1) > 0 )
                                    BF(x1==1) = wstar ./ (wstar + (1 - wstar) * exp(- diff(x1==1)));
                                end
                            end
                            
                            if ( isreal(BF) == 0 )
                                stop;
                            end
                            
                            sumgamma3{subj, 1}(chain1, :) = sumgamma3{subj, 1}(chain1, :) + BF;
                            count3{subj, 1}(chain1) = count3{subj, 1}(chain1) + 1;
                        end
                        
                        if ( logitzeta0(subj, chain1) < 0 )
                            zeta = (logita(subj) + logitb(subj) * exp(logitzeta0(subj, chain1))) ./ (1 + exp(logitzeta0(subj, chain1)));
                        else
                            zeta = (logita(subj) * exp(- logitzeta0(subj, chain1)) + logitb(subj)) ./ (exp(- logitzeta0(subj, chain1)) + 1);
                        end
                        probstar = sumgamma3{subj, 1}(chain1, :) / count3{subj, 1}(chain1);
                        
                        zetaAS = probstar .* zeta .* min(1./probstar, 1./(1 - probstar));
                        zetaDS = (1 - probstar) .* zeta .* min(1./probstar, 1./(1 - probstar));
                        
                        expec_change = sum(zetaAS .* (1 - probstar)) + sum(zetaDS .* probstar);
                        if ( expec_change < 1 )
                            logitzeta0(subj, chain1) = logitzeta0(subj, chain1) - log(expec_change);
                        end
                        
                        if ( logitzeta0(subj, chain1) < - 3 )
                            logitzeta0(subj, chain1) = - 3;
                        end
                        if ( logitzeta0(subj, chain1) > 3 )
                            logitzeta0(subj, chain1) = 3;
                        end
                        
                        if ( isreal(logitzeta0(subj, chain1)) == 0 )
                            expec_change
                            logitzeta0(subj, chain1)
                            stop;
                        end
                    end
                    
                    if ( size(C{subj, rep, chain}, 2) ~=  1 + sum(gamma{subj, 1}(rep, :, chain)))
                        stop;
                    end
    
                end
                
                
                
                
                if ( updateg == 1 )
                    newg = g(rep, chain) * exp(exp(loggsd(chain1)) * randn);
                    
                    logaccept = 0;
                    newloglike = zeros(1, numb_subjects);
                    newC = cell(numb_subjects, 1);
                    for subj = 1:numb_subjects
                        datastar = [ones(n(subj), 1) data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)];
                        if ( gprior == 1 )
                            n0star = zeros(sum(gamma(rep, :, chain)));
                            n0star(1:(sum(gamma(rep, :, chain)) + size(data_fixed, 2)), 2:(sum(gamma(rep, :, chain)) + size(data_fixed, 2))) = 1 / newg * datastar(:, 2:end)' * datastar(:, 2:end);
                        else
                            n0star = 1 / newg * eye(sum(gamma{subj, 1}(rep, :, chain)) + 1);
                            n0star(1, 1) = 1 / sigmasq_int(rep, chain);
                        end
                        newC{subj, 1} = inv(datastar' * datastar + n0star);
                        
                        XTy = [sum(resid{subj, 1}(:, rep, chain)); data{subj, 1}(:, gamma{subj, 1}(rep, :, chain)==1)' * resid{subj, 1}(:, rep, chain)];
                        newloglike(subj) = - 0.5 * sum(gamma{subj, 1}(rep, :, chain)) * log(newg);
                        newloglike(subj) = newloglike(subj) - 0.5 * log(det(datastar' * datastar + n0star));
                        newloglike(subj) = newloglike(subj) - (0.5 * n(subj) + sigmasqalpha(rep, chain)) * log(sigmasqbeta(rep, chain) + 0.5 * (resid{subj, 1}(:, rep, chain)' * resid{subj, 1}(:, rep, chain) - XTy' * newC{subj, 1} * XTy));
    
                        logaccept = logaccept + heat(chain1) * (newloglike(subj) - loglike(subj, rep, chain));
                    end
    
    %                 'here'
    %                 logaccept
                    logaccept = logaccept + 0.5 * log(newg) - log(1 + newg);
                    logaccept = logaccept - 0.5 * log(g(rep, chain)) + log(1 + g(rep, chain));
                    
                    accept = 1;
                    if ( isreal(logaccept) == 0 )
                        accept = 0;
                    elseif ( (isnan(logaccept) == 1) || (isinf(logaccept) == 1) )
                        accept = 0;
                    elseif ( logaccept < 0 )
                        accept = exp(logaccept);
                    end
                    
    %                 logaccept
    %                 accept
                    
                    if ( rand < accept )
                        g(rep, chain) = newg;
                        loglike(:, rep, chain) = newloglike;
                        for subj = 1:numb_subjects
                            C{subj, rep, chain} = newC{subj, 1};
                        end
                    end
                    
                    gaccept(chain1) = gaccept(chain1) + accept;
                    gcount(chain1) = gcount(chain1) + 1;
                    
                    loggsd(chain1) = loggsd(chain1) + 1 / it^0.55 * (accept - 0.234);
                end
                
            end
        end
        
        
        for rep = 1:numbofreps
            for chain = 1:numbofchains
                for subj = 1:numb_subjects
                    if ( size(C{subj, rep, chain}, 2) ~=  1 + sum(gamma{subj, 1}(rep, :, chain)))
                        rep
                        chain
                        subj
                        stop;
                    end
                end
            end
        end
        
        
        if ( length(heat) > 1 )
            
            % Swap move
            for rep = 1:numbofreps
                check = 1;
                while ( check == 1 )
                    chain1 = ceil(rand * length(heat));
                    check = (order(rep, chain1) == length(heat));
                end
                chain2 = find(order(rep, :) == order(rep, chain1) + 1);
                
                loglike1 = loglike(rep, chain1);
                loglike2 = loglike(rep, chain2);
                
                logaccept = heat(order(rep, chain1)) * loglike2 + heat(order(rep, chain2)) * loglike1 - heat(order(rep, chain1)) * loglike1 - heat(order(rep, chain2)) * loglike2;
                
                accept = 1;
                if ( logaccept  < 0 )
                    accept = exp(logaccept);
                end
                totalaccept2(order(rep, chain1)) = totalaccept2(order(rep, chain1)) + accept;
                totalcount2(order(rep, chain1)) = totalcount2(order(rep, chain1)) + 1;
                
                rhostar(order(rep, chain1)) = rhostar(order(rep, chain1)) + 1 / it^0.55 * (accept - 0.234);
                if ( rhostar(order(rep, chain1)) > 4 )
                    rhostar(order(rep, chain1)) = 4;
                end
                
                for j = 1:(length(heat)-1)
                    heat(j+1) = heat(j) * exp(- exp(rhostar(j)));
                end
                
                if ( rand < accept )
                    order(rep, chain1) = order(rep, chain1) + 1;
                    order(rep, chain2) = order(rep, chain2) - 1;
                end
            end
        end
        
        for rep = 1:numbofreps
            for chain = 1:numbofchains
                for subj = 1:numb_subjects
                    if ( size(C{subj, rep, chain}, 2) ~=  1 + sum(gamma{subj, 1}(rep, :, chain)))
                        stop;
                    end
                end
            end
        end
        
        if ( verLessThan('matlab', '9.7') == 0 )
            
            d.Value = it / (burnin+numbofits*thin);
            d.Message = 'Current Progress';
            
        else
            
            Val = floor(100 * it / (burnin+numbofits*thin));
            uialert(fig1, [num2str(Val) '% completed'], 'Fitting Run','Icon', 'info');
            
        end
    
        
        if ( (it > burnin) && (mod(it - burnin, thin) == 0) )
            
            for rep = 1:numbofreps
                top_chain = find(order(rep, :) == 1);
                it1 = (it - burnin) / thin;
                it2 = (it1 - 1) * numbofreps + rep;
                
                if ( it2 > length(holdsigmasqalpha) )
                    stop;
                end
                holdsigmasqalpha(it2) = sigmasqalpha(rep, top_chain);
                if ( it2 > length(holdsigmasqbeta) )
                    stop;
                end
                holdsigmasqbeta(it2) = sigmasqbeta(rep, top_chain);
                if ( it2 > length(holdg) )
                    stop;
                end
                holdg(it2) = g(rep, top_chain);
                holdg_fixed(it2) = g_fixed(rep, top_chain);
                
                for subj = 1:numb_subjects
                    
                    if ( size(C{subj, rep, top_chain}, 2) ~=  1 + sum(gamma{subj, 1}(rep, :, top_chain)))
                        stop;
                    end
                    
                    sumgamma{subj, 1} = sumgamma{subj, 1} + gamma{subj, 1}(rep, :, top_chain);
                    
    %                 if ( fixed == 1 )
    %                     holdloglike(subj, (it - burnin) / thin, rep) = loglike(subj, rep, top_chain) + sum(gamma{subj, 1}(rep, :, top_chain)) * log(w(rep, chain)) + sum(1-gamma{subj, 1}(rep, :, top_chain)) * log(1 - w(rep, chain));
    %                 else
    %                     holdloglike(subj, (it - burnin) / thin, rep) = loglike(subj, rep, top_chain) + gammaln(wa + sum(gamma{subj, 1}(rep, :, top_chain))) + gammaln(wb + p(subj) - sum(gamma{subj, 1}(rep, :, top_chain))) - gammaln(wa + wb + p(subj));
    %                 end
    %               holdaccept(subj, (it - burnin) / thin) = totalaccept(1) / totalcount(1);
    
                    if ( mode == 2 )
                        
                        datastar = [ones(n(subj), 1) data{subj, 1}(:, gamma{subj, 1}(rep, :, top_chain)==1)];
                        if ( gprior == 1 )
                            n0star = zeros(sum(gamma(rep, :, top_chain)));
                            n0star((1+size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma(rep, :, top_chain))), (1+size(data_fixed, 2)):(size(data_fixed, 2) + sum(gamma(rep, :, top_chain)))) = 1 / g(rep, top_chain) * datastar(:, 2:end)' * datastar(:, 2:end);
                        else
                            n0star = 1 / g(rep, top_chain) * eye(sum(gamma{subj, 1}(rep, :, top_chain)) + 1);
                            n0star(1, 1) = 1 / sigmasq_int(rep, chain);
                        end
                        
                        XTy = [sum(resid{subj, 1}(:, rep, top_chain)); data{subj, 1}(:, gamma{subj, 1}(rep, :, top_chain)==1)' * resid{subj, 1}(:, rep, top_chain)];
                        
                        alphastar = 0.5 * n(subj) + sigmasqalpha(rep, chain);
                        betastar = sigmasqbeta(rep, chain) + 0.5 * (resid{subj, 1}(:, rep, top_chain)' * resid{subj, 1}(:, rep, top_chain)) - 0.5 * XTy' * C{subj, rep, top_chain} * XTy;
    
    %                     size(XTy)
    %                     size(C{subj, rep, top_chain})
    %                     sum(gamma{subj, 1}(rep, :, top_chain))
    %                     size(n0star)
                        
                        sigmasq = 1 / gamrnd(alphastar, 1 / betastar);
                        
                        mustar = C{subj, rep, top_chain} * XTy;
                        varstar = sigmasq * C{subj, rep, top_chain};
                        
                        beta_local = mustar + chol(varstar)' * randn(size(C{subj, rep, top_chain}, 1), 1);
    
                        if ( subj > size(holdalpha, 1) )
                            stop;
                        end
                        if ( it2 > size(holdalpha, 2) )
                            stop;
                        end
                        holdalpha(subj, it2) = beta_local(1);
                        
                        if ( subj > size(holdgamma, 1) )
                            stop;
                        end
                        if ( it2 > size(holdgamma, 2) )
                            stop;
                        end
                        holdgamma{subj, it2}  = find(gamma{subj, 1}(rep, :, top_chain));
                        if ( subj > size(holdbeta, 1) )
                            stop;
                        end
                        if ( it2 > size(holdbeta, 2) )
                            stop;
                        end
                        holdbeta{subj, it2} =  beta_local(2:end);
    
                        if ( subj > size(holdgamma, 1) )
                            stop;
                        end
                        if ( it2 > size(holdgamma, 2) )
                            stop;
                        end
                        holdsigmasq(subj, it2) = sigmasq;
    
                    end
                end
                if ( it2 > size(holdbeta0, 1) )
                    stop;
                end
                holdbeta0(it2, :) = alpha_fixed(:, rep, top_chain);
                holdsigmasq_int(it2) = sigmasq_int(rep, top_chain);
            end
        end
    end
    
    for subj = 1:numb_subjects
        sumgamma{subj, 1} = sumgamma{subj, 1} / (numbofits * numbofreps);
        sumgamma2{subj, 1} = sumgamma2{subj, 1} / (numbofits * numbofreps);
    end
    
    zeta = zeros(numb_subjects, numbofchains);
    zetaAS = cell(numb_subjects, 1);
    zetaDS = cell(numb_subjects, 1);
    probstar = cell(numb_subjects, 1);
    for subj = 1:numb_subjects
        zetaAS{subj, 1} = zeros(numbofchains, p(subj));
        zetaDS{subj, 1} = zeros(numbofchains, p(subj));
        
        for chain1 = 1:numbofchains        
            if ( logitzeta0(subj, chain1) < 0 )
                zeta(subj, chain1) = (logita(subj) + logitb(subj) * exp(logitzeta0(subj, chain1))) ./ (1 + exp(logitzeta0(subj, chain1)));
            else
                zeta(subj, chain1) = (logita(subj) * exp(- logitzeta0(subj, chain1)) + logitb(subj)) ./ (exp(- logitzeta0(subj, chain1)) + 1);
            end
            
            probstar{subj, 1} = sumgamma3{subj, 1}(chain1, :) / count3{subj, 1}(chain1);
            
            zetaAS{subj, 1}(chain1, :) = probstar{subj, 1} .* zeta(subj, chain1) .* min(1./probstar{subj, 1}, 1./(1 - probstar{subj, 1}));
            zetaDS{subj, 1}(chain1, :) = (1 - probstar{subj, 1}) .* zeta(subj, chain1) .* min(1./probstar{subj, 1}, 1./(1 - probstar{subj, 1}));
        end
    end
    
    prob_inclusion = sumgamma;
    %modelsize = holdmodelsize;
    %gamma = holdgamma;
    alpha = holdalpha;
    beta = holdbeta;
    beta0 = holdbeta0;
    gamma = holdgamma;
    sigmasq = holdsigmasq;
    
    output = struct('sigmasq_int', holdsigmasq_int, 'g', holdg, 'sigmasqalpha', holdsigmasqalpha, 'sigmasqbeta', holdsigmasqbeta, 'zeta', zeta, 'zetaAS', zetaAS, 'zetaDS', zetaDS);