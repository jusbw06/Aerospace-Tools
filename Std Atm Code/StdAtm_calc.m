function [solution] = StdAtm_calc(Hp_vec,StdAtm) %General, allows
  %different models to be passed in
  T = zeros(1,length(Hp_vec));
  P = zeros(1,length(Hp_vec));
  rho_std = zeros(1,length(Hp_vec));
  a = zeros(1,length(Hp_vec));
  if length(StdAtm.layersHp)> length(StdAtm.A)
    ceiling = StdAtm.layersHp(length(StdAtm.layersHp));
  else
    ceiling = NaN;
  end
  for i = 1:length(Hp_vec)
    Hp = Hp_vec(i);
    ind = 1;
    for j = 1:length(StdAtm.layersHp)
      if Hp>StdAtm.layersHp(j)
        ind = j;
        if Hp> ceiling
          %disp('The entered height is more than the celiing for this model. Expect inaccuracy.')
          ind = j-1;
        end
      end
    end
    layer_handle = StdAtm.indcheck(ind);
    [T(i),P(i),rho_std(i),a(i)] = layer_handle(Hp,ind);
  end
  solution.T = T;
  solution.P = P;
  solution.rho_std = rho_std;
  solution.a = a;
end
