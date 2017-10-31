function R = ex_synth_swsat_tais(b)
% function R = ex_synth_swsat_tais(b)
% inputs: b (optative)
%
% outputs: bmax: max. feas. b
%
% Created: 07/fev/2010
% Updated: 11/fev/2010
% taisrc@dt.fee.unicamp.br

%Exemplos
exemplo=1;
fprintf('\n');

switch exemplo
    case 1
        %Exemplo 1
        if nargin <= 1
            [Ai,Bi,N,n,m] = example1(b);
            %R = synth_swsat_tais_d_yal(Ai,Bi);
            %R = synth_swsat_tais_finsler_d_yal(Ai,Bi);
            %R = synth_swsat_tais_slack_d_yal(Ai,Bi);
            %R = synth_swsat_MYG09_T1_d_yal(Ai,Bi);
            %rhoi = [5]; R = synth_swsat_tais_slack_otim_d_yal(Ai,Bi,rhoi);
            rhoi = 10; c = 0.2; r = 0.5;
            for rhoi = 1:9:20
                R = synth_swsat_tais_slack_poleplace_otim_d_yal(Ai,Bi,rhoi,c,r);
                if R.feas == 1
                    figure; hold on;
                    for i = 1:N
                        Wi = R.Wi(:,n*(i-1)+1:i*n);
                        Pi = inv(Wi);
                        projellisa2(Pi);
                        %V(x) = x'*P*x <= gamma => Vol=sqrt(det(inv(P*gamma)))                        
                        volume(i) = sqrt(det(Wi));
                    end
                    hold off;
                    title(['rhoi =',  num2str(rhoi), ',  volume =',  num2str(volume)]);
                    %display(rhoi);display(volume);disp('-----');
                    else
                        disp('Condição infactivel!');
                end
            end
        else
            disp('Exemplo 1 [Discrete swit. sys. w/ sat.]'); %Fonte: exemplo 1 [DH08]           
            disp('type');
            disp('-------------------')
            prec=0.1;            
            b=1.4-prec;
            while prec > 1e-4
                b=b+prec;
                [Ai,Bi] = example1(b);
                
                %R = synth_swsat_tais_d_yal(Ai,Bi);
                %R = synth_swsat_tais_finsler_d_yal(Ai,Bi);
                %R = synth_swsat_tais_slack_d_yal(Ai,Bi);
                %R = synth_swsat_MYG09_T1_d_yal(Ai,Bi);
                %rhoi = [10];
                %R = synth_swsat_tais_slack_otim_d_yal(Ai,Bi,rhoi);
                rhoi = [10]; c = 0.2; r = 0.7;
                R = synth_swsat_tais_slack_poleplace_otim_d_yal(Ai,Bi,rhoi,c,r);
                
                fprintf('b= %1.4f => %d [%2.2f],K=%d, L=%d\n',b,R.feas,R.cpusec,R.K,R.L)
                if R.feas == 0
                    b=b-prec;
                    prec=prec/10;
                end
            end
        end
        %Results: bmax
        %R = synth_swsat_tais_d_yal(Ai,Bi)          => b= 2.0650 [0.05],K=16, L=28
        %R = synth_swsat_tais_finsler_d_yal(Ai,Bi)  => b= 1.8770 [0.02],K=20, L=28
        %R = synth_swsat_tais_slack_d_yal(Ai,Bi)    => b= 2.0650 [0.03],K=24, L=28
        %R = synth_swsat_MYG09_T1_d_yal(Ai,Bi)      => b= 2.0650 [0.06],K=28, L=72
end