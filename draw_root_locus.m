function draw_root_locus(G)
% 打开本文件，在 MATLAB 编辑器中点“运行”，然后在命令窗口输入传递函数。
% 题目写 K* 时，输入 K/(s*(s+1)*(s+2))，其中 K 表示可变乘数。
% 也可输入 1/(s*(s+1)*(s+2))，程序自动扫描同一个乘数 K。
% 也可输入系数形式：tf([5 5 0], [1 1 1 -2])
% 需要 Control System Toolbox。适用于连续时间、实系数、无延迟 SISO 模型。

expression = '';
if nargin == 0
    s = tf('s'); %#ok<NASGU> % input 在当前函数中求值，输入时可直接使用 s。
    K = 1; % 输入表达式中的 K 仅作为乘法因子，取 1 得到根轨迹基准模型。
    expression = cleanPromptExpression(input('请输入传递函数（例如 K/(s*(s+1)*(s+2))，K 为可变乘数）：\nG = ','s'));
    try
        G = eval(expression);
    catch cause
        error('传递函数输入有误：%s\n可照这个格式输入：K/(s*(s+1)*(s+2))\nMATLAB 提示：%s', ...
            expression,cause.message);
    end
end
if ~isa(G, 'tf') || ~issiso(G) || G.Ts ~= 0 || hasdelay(G) || ~isproper(G)
    error('请输入连续时间、无延迟的单输入单输出有理传递函数。');
end
[N,D] = tfdata(G, 'v');
p = pole(G);
z = zero(G);
if isempty(p) || ~any(N) || ~isreal(N) || ~isreal(D)
    error('请输入实系数、非零且至少有一个极点的传递函数。');
end

% 首项系数之比保留在输入模型中。分子、分母各自首一化后，
% 1 + K*G(s) = 0 中的乘数对应 K* = gainScale * K。
gainScale = N(find(N~=0,1))/D(find(D~=0,1));
fprintf('根轨迹增益 K* = (%.8g) × K，K 为输入模型外乘的非负增益。\n',gainScale);
[w,K,continuousCase] = axisPoints(N,D);
% K=0 的位置由开环极点标记显示，不再重复标为虚轴交点。
positiveGain = K>0;
w = w(positiveGain);
K = K(positiveGain);
[breakS,breakK] = breakawayPoints(N,D,p,z);
[~,kg] = rlocus(G);
r = rlocus(G, unique([0;kg(:);K;breakK]).');

% 以零极点和虚轴交点为基准设置初始窗口，可用 MATLAB 工具栏缩放。
points = [0;p;z;1i*w;breakS];
q = numel(p)-numel(z);
if q > 0
    center = real((sum(p)-sum(z))/q);
    points = [points;center];
end
x = real(points);
y = imag(points);
xpad = max(1.4,0.28*(max(x)-min(x)));
ypad = max(1.5,0.23*(max(y)-min(y)));
xl = [min(x)-xpad,max(x)+1.35*xpad];
yl = [min(y)-ypad,max(y)+ypad];
step = 0.055*min(diff(xl),diff(yl));

fig = figure('Color','w','Name','根轨迹','NumberTitle','off', ...
    'WindowStyle','normal','Resize','on','MenuBar','figure','ToolBar','figure');
screen = get(groot,'ScreenSize');
fig.Position(3:4) = min(max(fig.Position(3:4),[930 580]),0.86*screen(3:4));
movegui(fig,'center');
top = showFormula(fig,N,D,'','',expression);
ax = axes(fig);
ax.Position = [0.09 0.12 0.82 top-0.19];
if ~isempty(breakS)
    ax.Position = [0.075 0.12 0.63 top-0.19];
end
hold(ax,'on');
axis(ax,'equal');
xlim(ax,xl);
ylim(ax,yl);
set(ax,'XAxisLocation','origin','YAxisLocation','origin','Box','off');
xlabel(ax,'实轴 Re(s)');
ylabel(ax,'虚轴 Im(s)');
title(ax,'根轨迹');

if q > 0
    radius = 2*max(hypot(max(abs(xl-center)),max(abs(yl))),1);
    angles = (2*(0:q-1)+(gainScale>0))*180/q;
    for angle = angles
        dx = cosd(angle);
        dy = sind(angle);
        plot(ax,center+[0 radius]*dx,[0 radius]*dy,'--', ...
            'Color',[0.65 0.65 0.65],'HandleVisibility','off');
    end
    fprintf('渐近线方向角（从正实轴逆时针量起）：');
    fprintf(' %.6g°',angles);
    fprintf('\n');
end

for k = 1:size(r,1)
    b = r(k,:);
    plot(ax,real(b),imag(b),'k-','LineWidth',1.8,'HandleVisibility','off');
    d = diff(b);
    inside = real(b)>xl(1)+step & real(b)<xl(2)-step & ...
        imag(b)>yl(1)+step & imag(b)<yl(2)-step;
    j = find(inside(1:end-1) & inside(2:end) & abs(d)>1e-5*step & abs(d)<2*step);
    if ~isempty(j)
        i = j(ceil(0.7*numel(j)));
        v = step*d(i)/abs(d(i));
        quiver(ax,real(b(i)),imag(b(i)),real(v),imag(v),0,'k', ...
            'LineWidth',1.4,'MaxHeadSize',1.5,'HandleVisibility','off');
    end
end

if q > 0
    drawAngleArcs(ax,center,angles,xl,yl);
end
markPoints(ax,p,'x',[0.75 0.1 0.1],-step,'开环极点');
markPoints(ax,z,'o',[0 0.35 0.7],step,'开环零点');
if ~isempty(breakS)
    color = [0.56 0.23 0];
    plot(ax,breakS,zeros(size(breakS)),'d','Color',color, ...
        'MarkerFaceColor',[1 0.74 0.35],'MarkerSize',7,'DisplayName','实轴分离点');
    panel = uipanel(fig,'Title','实轴分离点','FontSize',10, ...
        'BackgroundColor','w','Units','normalized','Position',[0.745 0.46 0.235 0.35]);
    rows = cell(numel(breakS),2);
    for k = 1:numel(breakS)
        rows(k,:) = {sprintf('(%.6g, 0)',breakS(k)),sprintf('%.6g',gainScale*breakK(k))};
    end
    uitable(panel,'Data',rows,'ColumnName',{'坐标 (x, y)','增益 K*'}, ...
        'RowName',[],'ColumnWidth',{92,58},'Units','normalized', ...
        'Position',[0.03 0.05 0.94 0.90]);
    disp('实轴分离点及其增益：');
    disp(gainTable(breakS,breakK,gainScale));
end
if ~isempty(w)
    plot(ax,zeros(size(w)),w,'ko', ...
        'MarkerFaceColor','k','MarkerSize',6,'DisplayName','虚轴交点');
    for k = 1:numel(w)
        text(ax,0.4*step,w(k)+0.65*step, ...
            sprintf('s = %+.4gj，K* = %.5g',w(k),gainScale*K(k)));
    end
    disp(gainTable(1i*w,K,gainScale));
elseif ~continuousCase
    disp('根轨迹与虚轴无 K>0 的有限增益交点。');
end
if continuousCase
    disp('此模型可能存在沿虚轴延伸的分支。当前仅单独检查原点，请结合图像查看。');
end
legend(ax,'show','Location','best');
end

function expression = cleanPromptExpression(raw)
% 接受从题目或命令窗口粘贴的常见数学写法。
expression = strtrim(raw);
from = {'（','）','［','］','＋','－','−','×','·','⋅','＊','÷','／','＾','＝','；','０','１','２','３','４','５','６','７','８','９'};
to = {'(',')','[',']','+','-','-','*','*','*','*','/','/','^','=',';','0','1','2','3','4','5','6','7','8','9'};
for k = 1:numel(from)
    expression = strrep(expression,from{k},to{k});
end
expression = regexprep(expression,'^\s*(?:>>\s*)?(?:[Gg](?:\s*\(\s*[sS]\s*\))?\s*=\s*)?','');
expression = regexprep(expression,'[;\\\s]+$','');
expression = strrep(expression,'**','^');
superscript = '⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺';
digits = '0123456789-+';
[first,last] = regexp(expression,'[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+');
for j = numel(first):-1:1
    [~,index] = ismember(expression(first(j):last(j)),superscript);
    exponent = digits(index);
    if startsWith(exponent,'-') || startsWith(exponent,'+')
        exponent = ['(' exponent ')'];
    end
    expression = [expression(1:first(j)-1) '^' exponent expression(last(j)+1:end)];
end
if ~isempty(regexp(expression,'\b(?:tf|zpk)\s*\(','once'))
    return; % 保留系数数组的括号和空格。
end
if contains(expression,'[') && ~isempty(regexp(expression,'\[[^\]]*[sS]','once'))
    expression = strrep(strrep(expression,'[','('),']',')');
end
expression = regexprep(expression,'\)\s*\(',')*(');
expression = regexprep(expression,'\)\s*(?=[sSKk])',')*');
expression = regexprep(expression,'(\d)\s*(?=[sSKk](?![A-Za-z0-9_]))','$1*');
expression = regexprep(expression,'(?<![A-Za-z0-9_.])(\d+(?:\.\d*)?)\s*\(','$1*(');
expression = regexprep(expression,'(?<![A-Za-z_])([sSKk])(?=\s*\()','$1*');
expression = regexprep(expression,'(?<![A-Za-z_])[Kk]\s*(?=[sS](?![A-Za-z0-9_]))','K*');
expression = regexprep(expression,'(?<![A-Za-z0-9_])s\s*(?=[Kk](?![A-Za-z0-9_]))','s*');
expression = regexprep(expression,'(?<![A-Za-z0-9_])S(?![A-Za-z0-9_])','s');
expression = regexprep(expression,'(?<![A-Za-z0-9_])k(?![A-Za-z0-9_])','K');
end

function top = showFormula(fig,numerator,denominator,prefix,suffix,expression)
% 表达式输入保留因子和括号。直接传 tf 对象时，从系数排版。
hasExpression = ~isempty(expression) && isempty(regexp(expression,'^\s*tf\s*\(','once'));
if hasExpression
    prefix = '';
    split = fractionSlash(expression);
    if split > 0
        num = expressionRows(expression(1:split-1));
        den = expressionRows(expression(split+1:end));
    else
        num = expressionRows(expression);
        den = {};
    end
else
    num = polynomialRows(numerator);
    den = polynomialRows(denominator);
end
rowCount = numel(num)+numel(den)+~isempty(suffix);
height = min(0.43,max(0.14,0.041*(rowCount+0.9)));
rowHeight = height/(rowCount+0.9);
top = 0.97-height;
fontSize = min(13,max(7,340*rowHeight));
annotation(fig,'rectangle',[0.065 top 0.89 height], ...
    'FaceColor',[0.97 0.98 0.99],'Color',[0.76 0.80 0.85], ...
    'Tag','TransferFunction');
x = 0.205;
width = 0.70;
barY = 0.97-(numel(num)+0.45)*rowHeight;
formulaText(fig,[0.075 barY-rowHeight/2 0.12 rowHeight], ...
    ['G(s) = ' prefix],fontSize);
for k = 1:numel(num)
    formulaText(fig,[x barY+(numel(num)-k)*rowHeight width rowHeight],num{k},fontSize);
end
if ~isempty(den)
    annotation(fig,'line',[x x+width],[barY barY],'Color',[0.23 0.27 0.31]);
end
for k = 1:numel(den)
    formulaText(fig,[x barY-k*rowHeight width rowHeight],den{k},fontSize);
end
if ~isempty(suffix)
    formulaText(fig,[x barY-(numel(den)+1)*rowHeight width rowHeight],suffix,fontSize);
end
end


function split = fractionSlash(expression)
% 只拆独立的最外层除法，以免改变混合运算的含义。
depth = 0;
split = 0;
for k = 1:numel(expression)
    ch = expression(k);
    if any(ch=='([{')
        depth = depth+1;
    elseif any(ch==')]}')
        depth = depth-1;
    elseif depth == 0
        previous = ' ';
        if k>1, previous = expression(k-1); end
        isExponentSign = any(previous=='eE') && k>2 && isstrprop(expression(k-2),'digit');
        if ch=='/'
            if split>0, split = 0; return; end
            split = k;
        elseif (any(ch=='+-') && ~isExponentSign && k>1 && ~any(previous=='^*/')) || (ch=='*' && split>0)
            split = 0;
            return;
        end
    end
end
end

function rows = expressionRows(expression)
% 在运算符前换行，再把输入的乘号和幂次换成数学字形。
line = regexprep(strtrim(expression),'\s+','');
rows = {};
while numel(line)>74
    positions = regexp(line(1:74),'[*+]','start');
    positions = positions(positions>15);
    if isempty(positions), break; end
    cut = positions(end);
    rows{end+1} = mathText(line(1:cut-1)); %#ok<AGROW>
    line = line(cut:end);
end
rows{end+1} = mathText(line);
end

function value = mathText(value)
value = strrep(value,'*',' \cdot ');
value = regexprep(value,'\^([+-]?\d+)','^{$1}');
end

function formulaText(fig,position,value,fontSize)
annotation(fig,'textbox',position,'String',value,'Interpreter','tex', ...
    'FitBoxToText','off','EdgeColor','none','Margin',0, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',fontSize,'Tag','TransferFunctionText');
end

function rows = polynomialRows(coefficients)
rows = {};
line = '';
for k = find(coefficients~=0)
    coefficient = coefficients(k);
    power = numel(coefficients)-k;
    factor = sprintf('%.8g',abs(coefficient));
    factor = regexprep(factor,'[eE]([+-]?)0*(\d+)','\\times 10^{$1$2}');
    if power > 0
        if abs(coefficient)==1, factor = ''; end
        factor = [factor 's']; %#ok<AGROW>
        if power>1, factor = [factor sprintf('^{%d}',power)]; end %#ok<AGROW>
    end
    if coefficient<0
        term = [' - ' factor];
    elseif isempty(rows) && isempty(line)
        term = factor;
    else
        term = [' + ' factor];
    end
    if ~isempty(line) && numel(line)+numel(term)>67
        rows{end+1} = strtrim(line); %#ok<AGROW>
        line = term;
    else
        line = [line term]; %#ok<AGROW>
    end
end
rows{end+1} = strtrim(line);
end

function data = gainTable(s,K,gainScale)
if gainScale==1
    data = table(s,K,'VariableNames',{'s','K_star'});
else
    data = table(s,K,gainScale*K,'VariableNames',{'s','K','K_star'});
end
end

function drawAngleArcs(ax,center,angles,xl,yl)
% 共轭渐近线共用上半平面的角度标注，实轴方向不标圆弧。
available = min([center-xl(1),xl(2)-center,-yl(1),yl(2)]);
arcColor = [0.32 0.32 0.32];
angles = angles(angles>0 & angles<180);
for angle = angles
    r = available*(0.20+0.32*angle/180);
    t = linspace(0,angle,50);
    plot(ax,center+r*cosd(t),r*sind(t),'-','Color',arcColor, ...
        'LineWidth',1,'HandleVisibility','off');
    mid = angle/2;
    labelRadius = r+0.10*available;
    text(ax,center+labelRadius*cosd(mid),labelRadius*sind(mid), ...
        sprintf('%.4g°',angle),'Color',arcColor, ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Clipping','on');
end
end

function [breakS,breakK] = breakawayPoints(N,D,p,z)
% 实轴分离点对应 K(s)=-D(s)/N(s) 的局部极大值。
% 分子 Q=D' N-D N' 在点的两侧由负变正时，K' 由正变负。
N = N(find(N~=0,1):end);
D = D(find(D~=0,1):end);
left = conv(polyder(D),N);
right = conv(D,polyder(N));
n = max(numel(left),numel(right));
Q = [zeros(1,n-numel(left)),left]-[zeros(1,n-numel(right)),right];
candidates = roots(Q);
realRoot = abs(imag(candidates))<=1e-7*max(1,abs(real(candidates)));
candidates = uniquetol(real(candidates(realRoot)),1e-8);
singularities = [p;z];
breakS = zeros(0,1);
breakK = zeros(0,1);
for x = candidates(:).'
    scale = max(1,abs(x));
    gap = min(abs(x-singularities));
    if gap <= 1e-7*scale
        continue;
    end
    gain = -polyval(D,x)/polyval(N,x);
    if ~isfinite(gain) || gain<=0
        continue;
    end
    other = candidates(candidates~=x);
    gap = min([gap;abs(other-x)]);
    dx = min(1e-3*scale,0.2*gap);
    if polyval(Q,x-dx)<0 && polyval(Q,x+dx)>0
        breakS(end+1,1) = x; %#ok<AGROW>
        breakK(end+1,1) = gain; %#ok<AGROW>
    end
end
end

function markPoints(ax,v,symbol,color,offset,label)
if isempty(v), return; end
plot(ax,real(v),imag(v),symbol,'Color',color,'MarkerSize',9, ...
    'LineWidth',1.8,'LineStyle','none','DisplayName',label);
for point = unique(v(:)).'
    text(ax,real(point)+0.25*abs(offset),imag(point)+offset, ...
        sprintf('(%.4g, %.4g)',real(point),imag(point)),'Color',color);
end
end

function [w,K,continuousCase] = axisPoints(N,D)
% D(jw)+K*N(jw)=0。使 D(jw)*conj(N(jw)) 的虚部为零，再回代求 K。
dn = (D/max(abs(D))).*(1i).^(numel(D)-1:-1:0);
nn = (N/max(abs(N))).*(-1i).^(numel(N)-1:-1:0);
product = imag(conv(dn,nn));
odd = product(end-1:-2:1);  % 提取 w 后得到关于 w^2 的多项式。
continuousCase = ~any(odd);
w = 0;
if ~continuousCase
    u = roots(fliplr(odd));
    positive = real(u)>0 & abs(imag(u))<=1e-8*abs(u);
    w = uniquetol([0;sqrt(real(u(positive)))],1e-8,'DataScale',1);
end
n = polyval(N,1i*w);
g = -polyval(D,1i*w)./n;
ok = abs(n)>10*eps*polyval(abs(N),w) & isfinite(g) & ...
    real(g)>=0 & abs(imag(g))<=1e-7*max(1,abs(real(g)));
w = w(ok);
K = real(g(ok));
K = [K(w>0);K];
w = [-w(w>0);w];
[w,order] = sort(w);
K = K(order);
end
