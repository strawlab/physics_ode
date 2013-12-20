load jd.mat
load jb.mat
load invI.mat
load invM.mat

nb=68;
damp_coef = 100.0;
dt = 0.001;
a = zeros(nb*6,nb*6);
for i=1:size(jb,1),
  b1_index = jb(i,1);
  b2_index = jb(i,2);

  a1_i_beg = (b1_index-1)*6;
  %a1_i_end = (b1_index)*6-1;
  a2_i_beg = (b2_index-1)*6;
  %a2_i_end = (b2_index)*6-1;

  %printf('filling out block %d %d\n',a1_i_beg,a2_i_beg);
  %fill in block
  for ai=1:6,
  for aj=1:6,
    a(a1_i_beg+ai,a1_i_beg+aj) = jd(i,ai)*damp_coef*jd(i,aj);
    a(a1_i_beg+ai,a2_i_beg+aj) = jd(i,ai)*damp_coef*jd(i,aj+6);
    a(a2_i_beg+ai,a1_i_beg+aj) = jd(i,aj+6)*damp_coef*jd(i,ai);
    a(a2_i_beg+ai,a2_i_beg+aj) = jd(i,aj+6)*damp_coef*jd(i,aj+6);
  end,
  end,

  %assemble inverse mass matrix
  m1 = zeros(6,6);
  m2 = zeros(6,6);
  m1(1,1) = invM(b1_index);
  m1(2,2) = invM(b1_index);
  m1(3,3) = invM(b1_index);
  m1(4:6,4:6) = invI( (b1_index-1)*3+1:(b1_index-1)*3+3, 1:3 );
  m2(1,1) = invM(b2_index);
  m2(2,2) = invM(b2_index);
  m2(3,3) = invM(b2_index);
  m2(4:6,4:6) = invI( (b2_index-1)*3+1:(b2_index-1)*3+3, 1:3 );

  %mass ratio of attached bodies
  m_check_12 = m1 * inv(m2);
  m_check_21 = m2 * inv(m1);
  m_cond(i) = max(abs(cond(m_check_12)),abs(cond(m_check_21))) ;

  %m1
  %a(a1_i_beg:a1_i_beg+5,a1_i_beg:a1_i_beg+5)
  %multiply by invM
  a(a1_i_beg+1:a1_i_beg+6,a1_i_beg+1:a1_i_beg+6) = dt*m1*a(a1_i_beg+1:a1_i_beg+6,a1_i_beg+1:a1_i_beg+6);
  a(a1_i_beg+1:a1_i_beg+6,a2_i_beg+1:a2_i_beg+6) = dt*m1*a(a1_i_beg+1:a1_i_beg+6,a2_i_beg+1:a2_i_beg+6);
  a(a2_i_beg+1:a2_i_beg+6,a1_i_beg+1:a1_i_beg+6) = dt*m2*a(a2_i_beg+1:a2_i_beg+6,a1_i_beg+1:a1_i_beg+6);
  a(a2_i_beg+1:a2_i_beg+6,a2_i_beg+1:a2_i_beg+6) = dt*m2*a(a2_i_beg+1:a2_i_beg+6,a2_i_beg+1:a2_i_beg+6);

end,
contour(a)

ident = eye(nb*6,nb*6);
a_expansion = ident + a - a*a + a*a*a - a*a*a*a;

E = (ident + a);
Einv = inv(E);

contour(Einv);
