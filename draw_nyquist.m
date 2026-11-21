function draw_nyquist(G,showNegative)
%DRAW_NYQUIST 输入传递函数，生成正频率 Nyquist 图。
% 在命令窗口运行 draw_nyquist，然后输入式子。含 K 时，再输入 K 的数值。
% draw_nyquist(G,true) 同时显示负频率。需要 Control System Toolbox。
expression = '';
gainCaption = '';
if nargin<2, showNegative = false; end
if nargin==0
    s = tf('s'); %#ok<NASGU>
    expression = cleanPromptExpression(input('请输入传递函数（例如 K/(s*(s+1)*(s+2))）：\nG = ','s'));
    if ~isempty(regexp(expression,'(?<!\w)K(?!\w)','once'))
        K = input('请输入增益 K 的数值：\nK = ');
        validateattributes(K,{'numeric'},{'real','scalar','finite','nonzero'},mfilename,'K');
        gainCaption = sprintf('K = %.8g',K);
    end
    try
        G = eval(expression);
    catch cause
        error('传递函数输入有误：%s\n可照这个格式输入：K/(s*(s+1)*(s+2))\nMATLAB 提示：%s', ...
            expression,cause.message);
    end
end
if ~isa(G,'tf') || ~issiso(G) || G.Ts~=0 || ~isreal(G) || ~isproper(G)
    error('请输入连续时间、实系数、正则有理的单输入单输出传递函数。');
end
[num,den] = tfdata(G,'v');
if ~any(num), error('请输入非零传递函数。'); end

m = allmargin(G);
[wc,pm] = finiteCrossings(m.PMFrequency,m.PhaseMargin);
[wg,gm] = finiteCrossings(m.GMFrequency,m.GainMargin);
% 单位圆交点在候选频率附近用实际响应细化。
for k = 1:numel(wc)
    level = @(f) abs(evalfr(G,1i*f))-1;
    bracket = wc(k)*[0.999 1.001];
    values = arrayfun(level,bracket);
    if all(isfinite(values)) && sign(values(1))*sign(values(2))<0
        wc(k) = fzero(level,bracket);
    end
end
gc = responseAt(G,wc);
gg = responseAt(G,wg);
[~,~,autoW] = bode(G);
autoW = autoW(isfinite(autoW) & autoW>0);
lo = min([0.01;autoW(:);wc/10;wg/10])/10;
hi = max([100;autoW(:);wc*10;wg*10])*10;
w = unique([logspace(log10(lo),log10(hi),3600).';wc;wg]);
r = responseAt(G,w);

% 直流端与高频极限由模型确定。原点有净极点时，低频端趋于无穷。
n0 = find(num~=0,1,'last');
d0 = find(den~=0,1,'last');
originPower = (numel(num)-n0)-(numel(den)-d0);
dc = NaN;
if originPower>0, dc = 0; end
if originPower==0, dc = num(n0)/den(d0); end
n1 = find(num~=0,1);
d1 = find(den~=0,1);
delay = G.InputDelay+G.OutputDelay+G.IODelay;
hf = NaN;
if numel(num)-n1<numel(den)-d1
    hf = 0;
elseif delay==0
    hf = num(n1)/den(d1);
end

fig = figure('Color','w','Name','Nyquist 图','NumberTitle','off', ...
    'WindowStyle','normal','Resize','on','MenuBar','figure','ToolBar','figure');
screen = get(groot,'ScreenSize');
fig.Position(3:4) = min(max(fig.Position(3:4),[930 580]),0.86*screen(3:4));
movegui(fig,'center');
prefix = '';
if delay>0, prefix = sprintf('e^{-(%.8g)s}',delay); end
top = showFormula(fig,num,den,prefix,gainCaption,expression);
ax = axes(fig,'Position',[0.11 0.12 0.78 top-0.19]);
hold(ax,'on');
axis(ax,'equal');
key = [1;abs(dc);abs(hf);abs(gc);abs(gg)];
radius = max(1.6,1.25*max(key(isfinite(key))));
xlim(ax,[-radius radius]); ylim(ax,[-radius radius]);
set(ax,'XAxisLocation','origin','YAxisLocation','origin','Box','off');
xlabel(ax,'实轴 Re\{G(j\omega)\}');
ylabel(ax,'虚轴 Im\{G(j\omega)\}');
title(ax,'Nyquist 图');
grid(ax,'on');

t = linspace(0,2*pi,400);
plot(ax,cos(t),sin(t),'--','Color',[0.55 0.55 0.55], ...
    'DisplayName','单位圆 |G| = 1');
if showNegative
    plot(ax,real(r),-imag(r),'--','Color',[0.66 0.72 0.80], ...
        'DisplayName','负频率');
end
blue = [0 0.29 0.64];
plot(ax,real(r),imag(r),'Color',blue,'LineWidth',1.9, ...
    'DisplayName','正频率（箭头指向频率增大）');
delta = diff(r);
visible = max(abs(real(r)),abs(imag(r)))<0.86*radius & abs(r)>0.05*radius;
ids = find(visible(1:end-1) & abs(delta)>1e-5*radius & abs(delta)<0.08*radius);
if ~isempty(ids)
    j = ids(ceil(numel(ids)/2));
    arrow = 0.115*radius*delta(j)/abs(delta(j));
    quiver(ax,real(r(j)),imag(r(j)),real(arrow),imag(arrow),0, ...
        'Color',blue,'LineWidth',1.5,'MaxHeadSize',1.4,'HandleVisibility','off');
end
plot(ax,-1,0,'rx','MarkerSize',10,'LineWidth',2,'DisplayName','临界点 (-1,0)');
if isfinite(dc) && dc==hf
    markEnd(ax,dc,'\omega = 0, \infty');
else
    markEnd(ax,dc,'\omega = 0');
    markEnd(ax,hf,'\omega \rightarrow \infty');
end
if ~isfinite(dc)
    text(ax,0.02,0.03,'低频端 \omega \rightarrow 0：|G| \rightarrow \infty', ...
        'Units','normalized','Color',blue);
end

orange = [0.74 0.25 0.06];
for k = 1:numel(wc)
    plot(ax,[0 real(gc(k))],[0 imag(gc(k))],'-','Color',orange,'HandleVisibility','off');
    plot(ax,real(gc(k)),imag(gc(k)),'o','Color',orange,'MarkerFaceColor',orange, ...
        'HandleVisibility','off');
    text(ax,real(gc(k)),0.035*radius,sprintf('  wc = %.3f',wc(k)), ...
        'Color',orange,'VerticalAlignment','bottom','Tag','NyquistWc');
    % 以负实轴为参考，取通往该单位圆交点的短弧。
    offset = mod(angle(gc(k))+2*pi,2*pi)-pi;
    a = linspace(-pi,-pi+offset,60);
    arcRadius = min(0.61,0.36+0.07*(k-1));
    plot(ax,arcRadius*cos(a),arcRadius*sin(a),'Color',orange,'HandleVisibility','off');
    mid = -pi+offset/2;
    text(ax,(arcRadius+0.08)*cos(mid),(arcRadius+0.08)*sin(mid), ...
        sprintf('\\gamma = %.4g^{\\circ}',pm(k)),'Color',orange, ...
        'HorizontalAlignment','center','Tag','NyquistPhaseMargin');
end
green = [0 0.46 0.28];
for k = 1:numel(wg)
    plot(ax,real(gg(k)),0,'d','Color',green,'MarkerFaceColor',green,'HandleVisibility','off');
    text(ax,real(gg(k)),0,sprintf('  wg = %.3f',wg(k)), ...
        'Color',green,'VerticalAlignment','top','Tag','NyquistWg');
end
legend(ax,'show','Location','best');
if isempty(wc), disp('没有有限正频率的单位圆交点 wc。'); end
if isempty(wg), disp('没有有限正频率的负实轴交点 wg。'); end
if ~isempty(wc)
    disp(table(wc,real(gc),imag(gc),pm,'VariableNames',{'wc_rad_s','Real','Imag','PhaseMargin_deg'}));
end
if ~isempty(wg)
    disp(table(wg,real(gg),abs(gg),gm,20*log10(gm), ...
        'VariableNames',{'wg_rad_s','Real','InverseGainMargin','GainMargin','GainMargin_dB'}));
end
end

function [frequency,margin] = finiteCrossings(frequency,margin)
frequency = frequency(:);
margin = margin(:);
keep = frequency>0 & isfinite(frequency) & isfinite(margin);
frequency = frequency(keep);
margin = margin(keep);
end

function response = responseAt(G,w)
if isempty(w)
    response = zeros(0,1);
else
    response = reshape(freqresp(G,w),[],1);
end
end

function markEnd(ax,point,label)
if isfinite(point)
    plot(ax,point,0,'ko','MarkerSize',4,'HandleVisibility','off');
    text(ax,point,0,['  ' label],'VerticalAlignment','bottom');
end
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
