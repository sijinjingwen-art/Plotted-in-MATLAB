function draw_bode(G)
%DRAW_BODE 在同一张图上画幅值和相角，并标出 0 dB 交越点 wc。
% 运行本文件，在命令窗口 G = 后输入表达式，例如 1/(s*(s+1))。
% 表达式含可变增益 K 时，例如 K/(s*(s+1))，接着输入 K 的数值。
% 需要 Control System Toolbox，输入连续时间单输入单输出传递函数。

expression = '';
gainCaption = '';
if nargin == 0
    s = tf('s'); %#ok<NASGU> % 供传递函数表达式使用。
    expression = cleanPromptExpression(input('请输入传递函数（例如 K/(s*(s+1))）：\nG = ','s'));
    if ~isempty(regexp(expression,'(?<!\w)K(?!\w)','once'))
        K = input('请输入增益 K 的数值：\nK = ');
        validateattributes(K,{'numeric'},{'real','scalar','finite','nonzero'}, ...
            mfilename,'K');
        gainCaption = sprintf('K = %.8g',K);
    end
    try
        G = eval(expression);
    catch cause
        error('传递函数输入有误：%s\n可照这个格式输入：K/(s*(s+1)*(s+2))\nMATLAB 提示：%s', ...
            expression,cause.message);
    end
end
if ~isa(G,'tf') || ~issiso(G) || G.Ts ~= 0 || ~isreal(G)
    error('请输入连续时间、实系数的单输入单输出传递函数。');
end
[num,den] = tfdata(G,'v');
if ~any(num)
    error('请输入非零传递函数。零系统的相角没有定义。');
end

% 检查无延迟模型是否在所有频率上均为 0 dB，例如全通系统。
mirrorNum = num.*(-1).^(numel(num)-1:-1:0);
mirrorDen = den.*(-1).^(numel(den)-1:-1:0);
a = conv(num,mirrorNum);
b = conv(den,mirrorDen);
n = max(numel(a),numel(b));
a = [zeros(1,n-numel(a)),a];
b = [zeros(1,n-numel(b)),b];
allZeroDb = isequal(a,b);
wc = zeros(0,1);
if ~allZeroDb
    margins = allmargin(G);
    candidates = margins.PMFrequency(:);
    wc = unique(candidates(isfinite(candidates) & candidates>0));
    % 在候选点附近用实际频率响应细化穿越频率。
    level = @(f) abs(evalfr(G,1i*f))-1;
    for k = 1:numel(wc)
        bracket = wc(k)*[0.999 1.001];
        levels = arrayfun(level,bracket);
        if all(isfinite(levels)) && sign(levels(1))*sign(levels(2))<0
            wc(k) = fzero(level,bracket);
        end
    end
    wc = unique(wc);
end

[~,~,autoW] = bode(G);
[v,lowGainDb,corners,slopeChange] = asymptoteData(G,num,den);
autoW = autoW(isfinite(autoW) & autoW>0);
if isempty(autoW)
    autoW = [0.01;100];
end
lo = min([autoW(:);wc/2;corners/2]);
hi = max([autoW(:);wc*2;corners*2]);
if lo == hi
    lo = lo/10;
    hi = hi*10;
end
w = unique([logspace(log10(lo),log10(hi),2000).';wc;corners]);
[mag,phase] = bode(G,w);
magDb = 20*log10(mag(:));
phase = phase(:);
[~,idx] = ismember(wc,w);
phaseAtWc = phase(idx);
edges = [lo;corners;hi];
asymDb = asymptoteDb(edges,v,lowGainDb,corners,slopeChange);
slopes = -20*v+[0;cumsum(slopeChange)];

% 按 20 dB 对应 90° 对齐两根纵轴，并扩展范围容纳两条曲线。
phasePerDb = 90/20;
finiteMag = [magDb(isfinite(magDb));asymDb];
finitePhase = phase(isfinite(phase));
leftLow = 20*floor(min([finiteMag;finitePhase/phasePerDb;-40])/20);
leftHigh = 20*ceil(max([finiteMag;finitePhase/phasePerDb;40])/20);
leftTicks = leftLow:20:leftHigh;

fig = figure('Color','w','Name','Bode 图','NumberTitle','off', ...
    'WindowStyle','normal','Resize','on','MenuBar','figure','ToolBar','figure');
screen = get(groot,'ScreenSize');
fig.Position(3:4) = min(max(fig.Position(3:4),[930 580]),0.86*screen(3:4));
movegui(fig,'center');
delay = G.InputDelay+G.OutputDelay+G.IODelay;
delayPrefix = '';
if delay > 0
    delayPrefix = sprintf('e^{-(%.8g)s}',delay);
end
top = showFormula(fig,num,den,delayPrefix,gainCaption,expression);
ax = axes(fig);
ax.Position = [0.10 0.12 0.80 top-0.19];
ax.XScale = 'log';
hold(ax,'on');
blue = [0 0.34 0.70];
orange = [0.78 0.25 0.08];
green = [0 0.44 0.31];
yyaxis(ax,'left');
gainLine = plot(ax,w,magDb,'Color',blue,'LineWidth',1.7);
asymLine = plot(ax,edges,asymDb,'-.','Color',green,'LineWidth',1.7, ...
    'Tag','MagnitudeAsymptote');
ylabel(ax,'模值 20 log_{10}|G(j\omega)| (dB)');
yline(ax,0,'k-','LineWidth',2.2,'HandleVisibility','off');
ylim(ax,[leftLow leftHigh]);
yticks(ax,leftTicks);
for k = 1:numel(corners)
    xline(ax,corners(k),'k--','LineWidth',1.1, ...
        'HandleVisibility','off','Tag','BodeCorner');
    plot(ax,corners(k),asymDb(k+1),'d','Color',green, ...
        'MarkerFaceColor','w','HandleVisibility','off');
    if mod(k,2)==1, side = 'top'; else, side = 'bottom'; end
    text(ax,corners(k),0,sprintf('w = %.5g',corners(k)), ...
        'Color',green,'HorizontalAlignment','center', ...
        'VerticalAlignment',side,'Tag','BodeCornerLabel');
end
for k = 1:numel(slopes)
    mid = exp((log(edges(k))+log(edges(k+1)))/2);
    slopeLabel = sprintf('%+g dB/dec',slopes(k));
    if slopes(k)==0, slopeLabel = '0 dB/dec'; end
    labelY = slopeLabelHeight(mid,mean(asymDb(k:k+1)),w,magDb, ...
        phase/phasePerDb,edges,asymDb,[leftLow leftHigh]);
    text(ax,mid,labelY,slopeLabel,'Color',green, ...
        'HorizontalAlignment','center','VerticalAlignment','middle','Tag','BodeSlope');
end
for k = 1:numel(wc)
    xline(ax,wc(k),':','Color',[0.55 0.55 0.55],'HandleVisibility','off');
    plot(ax,wc(k),0,'o','Color',blue,'MarkerFaceColor',blue,'HandleVisibility','off');
    text(ax,wc(k),0,sprintf('  wc = %.5g',wc(k)), ...
        'Color',blue,'VerticalAlignment','bottom');
end

yyaxis(ax,'right');
phaseLine = plot(ax,w,phase,'--','Color',orange,'LineWidth',1.7);
ylabel(ax,'相角 ∠G(j\omega) (°)');
ylim(ax,phasePerDb*[leftLow leftHigh]);
yticks(ax,phasePerDb*leftTicks);
for k = 1:numel(wc)
    plot(ax,wc(k),phaseAtWc(k),'s','Color',orange, ...
        'MarkerFaceColor',orange,'HandleVisibility','off');
    text(ax,wc(k),phaseAtWc(k),sprintf('  (%.5g, %.4g°)',wc(k),phaseAtWc(k)), ...
        'Color',orange,'VerticalAlignment','top');
end
ax.YAxis(1).Color = blue;
ax.YAxis(2).Color = orange;
xlabel(ax,'角频率 \omega (rad/s)');
title(ax,'Bode 图');
grid(ax,'on');
xlim(ax,[lo hi]);
legend(ax,[gainLine asymLine phaseLine],{'模值','幅值渐近线','相角'},'Location','best');

fprintf('原点净极点数 v = %d，低频渐近线斜率 = %g dB/dec。\n',v,slopes(1));
if ~isempty(corners)
    disp(table(corners,asymDb(2:end-1),slopes(2:end), ...
        'VariableNames',{'Corner_rad_s','Asymptote_dB','SlopeAfter_dB_dec'}));
end
if allZeroDb
    disp('模值始终为 0 dB，wc 无唯一值。');
elseif isempty(wc)
    disp('没有正频率的 0 dB 交越点。');
else
    disp(table(wc,zeros(size(wc)),phaseAtWc, ...
        'VariableNames',{'wc_rad_s','Magnitude_dB','Phase_deg'}));
end
end

function y = slopeLabelHeight(wMid,base,w,magDb,phaseOnLeft,edges,asymDb,limits)
% 在标签左右取样，根据到幅值、相角、渐近线和 0 dB 线的距离选位置。
span = diff(limits);
logW = log10(w);
halfWidth = 0.038*diff(log10(edges([1 end])));
probes = log10(wMid)+linspace(-halfWidth,halfWidth,9);
probes = min(max(probes,logW(1)),logW(end));
nearMag = interp1(logW,magDb,probes);
nearPhase = interp1(logW,phaseOnLeft,probes);
nearAsym = interp1(log10(edges),asymDb,probes);
nearCurves = [0 nearMag nearPhase nearAsym];
nearCurves = nearCurves(isfinite(nearCurves));
candidates = base+span*[0.055 -0.055 0.09 -0.09 0.13 -0.13];
candidates = min(max(candidates,limits(1)+0.028*span),limits(2)-0.028*span);
dist = min(abs(candidates(:)-nearCurves),[],2);
score = min(dist,0.065*span)-0.22*abs(candidates(:)-base);
[~,pick] = max(score);
y = candidates(pick);
end

function [v,lowGainDb,corners,change] = asymptoteData(G,num,den)
% v 是原点极点数减去原点零点数。低频系数取去掉 s 因子后的常数项之比。
iz = find(num~=0,1,'last');
ip = find(den~=0,1,'last');
v = (numel(den)-ip)-(numel(num)-iz);
lowGainDb = 20*(log10(abs(num(iz)))-log10(abs(den(ip))));
z = abs(zero(G));
p = abs(pole(G));
freq = [z(:);p(:)];
contribution = [ones(numel(z),1);-ones(numel(p),1)];
valid = freq>0 & isfinite(freq);
if ~any(valid)
    corners = zeros(0,1);
    change = zeros(0,1);
    return;
end
% 对数频率上合并同一转折点的零极点，复共轭对自然计两次。
[logCorners,~,group] = uniquetol(log10(freq(valid)),1e-6,'DataScale',1);
change = 20*accumarray(group,contribution(valid));
keep = change~=0;
corners = 10.^logCorners(keep);
change = change(keep);
end

function db = asymptoteDb(w,v,lowGainDb,corners,change)
logW = log10(w);
db = lowGainDb-20*v*logW;
if ~isempty(corners)
    db = db+max(0,logW-log10(corners).')*change;
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
