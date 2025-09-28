import java.awt.*;
import java.awt.geom.Path2D;
import java.io.*;
import java.nio.file.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.List;
import javax.swing.*;

public class App {
    public interface Function1D {
        double f(double x, EvalContext ctx);
        default double xmin(EvalContext ctx) { return Double.NEGATIVE_INFINITY; }
        default double xmax(EvalContext ctx) { return Double.POSITIVE_INFINITY; }
        default Function1D symbolicDerivative() { return null; }
        default String pretty() { return getClass().getSimpleName(); }
    }

    public static final class EvalContext {
        public final Map<String, Double> params;
        public EvalContext(Map<String, Double> params) { this.params = params; }
        public double p(String name, double def) {
            Double v = params.get(name);
            return v != null ? v : def;
        }
    }

    public interface DifferentiationScheme {
        double df(Function1D f, double x, double h, EvalContext ctx);
    }

    public static final class FivePointCentral implements DifferentiationScheme {
        private static final double EPS = 1e-8;
        @Override public double df(Function1D f, double x, double h, EvalContext ctx) {
            double xmin = f.xmin(ctx), xmax = f.xmax(ctx);
            double hh = adaptH(x, h, xmin, xmax);
            boolean can5 = (x - 2*hh >= xmin || Double.isInfinite(xmin)) && (x + 2*hh <= xmax || Double.isInfinite(xmax));
            if (can5) {
                double fm2 = f.f(x - 2*hh, ctx);
                double fm1 = f.f(x - hh, ctx);
                double fp1 = f.f(x + hh, ctx);
                double fp2 = f.f(x + 2*hh, ctx);
                return (-fp2 + 8*fp1 - 8*fm1 + fm2) / (12*hh);
            }
            boolean canF = (x + 2*hh <= xmax || Double.isInfinite(xmax));
            boolean canB = (x - 2*hh >= xmin || Double.isInfinite(xmin));
            if (canF) {
                double f0 = f.f(x, ctx);
                double fp1 = f.f(x + hh, ctx);
                double fp2 = f.f(x + 2*hh, ctx);
                return (-3*f0 + 4*fp1 - fp2) / (2*hh);
            }
            if (canB) {
                double f0 = f.f(x, ctx);
                double fm1 = f.f(x - hh, ctx);
                double fm2 = f.f(x - 2*hh, ctx);
                return (3*f0 - 4*fm1 + fm2) / (2*hh);
            }
            double tiny = Math.max(EPS, hh * 0.5);
            return (f.f(x + tiny, ctx) - f.f(x - tiny, ctx)) / (2*tiny);
        }
        private static double adaptH(double x, double h, double xmin, double xmax) {
            double hh = h;
            if (!Double.isInfinite(xmin) || !Double.isInfinite(xmax)) {
                double left = x - xmin;
                double right = xmax - x;
                double hMax = Math.min(
                        Double.isInfinite(left) ? Double.POSITIVE_INFINITY : left/2.1,
                        Double.isInfinite(right)? Double.POSITIVE_INFINITY : right/2.1
                );
                hh = Math.min(h, Math.max(1e-8, hMax));
            }
            return hh;
        }
    }

    public static abstract class Expr {
        abstract double eval(double x, EvalContext ctx);
        abstract Expr d();
        abstract String pretty();
        static Expr C(double v) { return new Const(v); }
        static Expr X() { return new Var("x"); }
        static Expr P(String name) { return new Param(name); }
        static Expr add(Expr a, Expr b) { return new Bin("+", a, b); }
        static Expr sub(Expr a, Expr b) { return new Bin("-", a, b); }
        static Expr mul(Expr a, Expr b) { return new Bin("*", a, b); }
        static Expr div(Expr a, Expr b) { return new Bin("/", a, b); }
        static Expr pow(Expr a, Expr b) { return new Bin("^", a, b); }
        static Expr sin(Expr a) { return new Uni("sin", a); }
        static Expr cos(Expr a) { return new Uni("cos", a); }
        static Expr exp(Expr a) { return new Uni("exp", a); }
        static Expr ln(Expr a) { return new Ln(a); }
    }

    public static final class Const extends Expr {
        final double v; Const(double v){this.v=v;}
        double eval(double x, EvalContext ctx){return v;}
        Expr d(){return C(0);}
        String pretty(){return format(v);}
    }

    public static final class Var extends Expr {
        final String name; Var(String n){this.name=n;}
        double eval(double x, EvalContext ctx){return name.equals("x")?x:ctx.p(name,0);}
        Expr d(){return name.equals("x")?C(1):C(0);}
        String pretty(){return name;}
    }

    public static final class Param extends Expr {
        final String name; Param(String n){this.name=n;}
        double eval(double x, EvalContext ctx){return ctx.p(name,0);}
        Expr d(){return C(0);}
        String pretty(){return name;}
    }

    public static final class Bin extends Expr {
        final String op; final Expr a,b;
        Bin(String op, Expr a, Expr b){this.op=op; this.a=a; this.b=b;}
        double eval(double x, EvalContext ctx){
            double u=a.eval(x,ctx), v=b.eval(x,ctx);
            return switch(op){
                case "+"->u+v; case "-"->u-v; case "*"->u*v; case "/"->u/v; case "^"->Math.pow(u,v);
                default->Double.NaN;
            };
        }
        Expr d(){
            return switch(op){
                case "+" -> add(a.d(), b.d());
                case "-" -> sub(a.d(), b.d());
                case "*" -> add(mul(a.d(), b), mul(a, b.d()));
                case "/" -> div(sub(mul(a.d(), b), mul(a, b.d())), pow(b, C(2)));
                case "^" -> {
                    Expr u=a, v=b;
                    yield mul(pow(u,v), add(
                            mul(v.d(), new Ln(u)),
                            mul(v, div(u.d(), u))
                    ));
                }
                default -> C(Double.NaN);
            };
        }
        String pretty(){return "("+a.pretty()+" "+op+" "+b.pretty()+")";}
    }

    public static final class Uni extends Expr {
        final String op; final Expr a;
        Uni(String op, Expr a){this.op=op; this.a=a;}
        double eval(double x, EvalContext ctx){
            double u=a.eval(x,ctx);
            return switch(op){
                case "sin"->Math.sin(u);
                case "cos"->Math.cos(u);
                case "exp"->Math.exp(u);
                default->Double.NaN;
            };
        }
        Expr d(){
            return switch(op){
                case "sin" -> mul(cos(a), a.d());
                case "cos" -> mul(C(-1), mul(sin(a), a.d()));
                case "exp" -> mul(this, a.d());
                default -> C(Double.NaN);
            };
        }
        String pretty(){return op+"("+a.pretty()+")";}
    }

    public static final class Ln extends Expr {
        final Expr a; Ln(Expr a){this.a=a;}
        double eval(double x, EvalContext ctx){ return Math.log(a.eval(x,ctx)); }
        Expr d(){ return div(a.d(), a); }
        String pretty(){ return "ln("+a.pretty()+")"; }
    }

    // ----------------- Parser & AnalyticFunction -----------------
    public static final class Parser {
        private final String s; private int i=0;
        public Parser(String s){ this.s=s.replaceAll("\\s+", ""); }
        public Expr parse(){
            List<Token> out = new ArrayList<>();
            Deque<Token> ops = new ArrayDeque<>();
            Token prev = null;
            while (i < s.length()){
                char c = s.charAt(i);
                if (Character.isDigit(c) || c=='.') { out.add(readNumber()); prev = out.get(out.size()-1); continue; }
                if (Character.isLetter(c)) { Token t = readIdentOrFunc();
                    if (t.type==T.FUNC) ops.push(t); else out.add(t);
                    prev = t; continue; }
                if (c=='('){ ops.push(new Token(T.LP,"(")); i++; prev=null; continue; }
                if (c==')'){
                    i++; while(!ops.isEmpty() && ops.peek().type!=T.LP) out.add(ops.pop());
                    if (ops.isEmpty()) throw new RuntimeException("Mismatched )");
                    ops.pop();
                    if (!ops.isEmpty() && ops.peek().type==T.FUNC) out.add(ops.pop());
                    prev=null; continue; }
                if ("+-*/^".indexOf(c)>=0){
                    String op = String.valueOf(c); i++;
                    if (op.equals("-") && (prev==null || prev.type==T.OP || prev.type==T.LP)){
                        out.add(new Token(T.NUM, "0"));
                        Token t = new Token(T.OP, "-");
                        while(!ops.isEmpty() && ops.peek().type==T.OP && prec(ops.peek().text) >= prec(t.text)) out.add(ops.pop());
                        ops.push(t); prev=t; continue;
                    }
                    Token t = new Token(T.OP, op);
                    while(!ops.isEmpty() && ops.peek().type==T.OP &&
                            ((rightAssoc(op)) ? (prec(ops.peek().text) > prec(op)) : (prec(ops.peek().text) >= prec(op)))) out.add(ops.pop());
                    ops.push(t); prev=t; continue;
                }
                throw new RuntimeException("Unexpected char: "+c);
            }
            while(!ops.isEmpty()){
                Token t = ops.pop();
                if (t.type==T.LP) throw new RuntimeException("Mismatched (");
                out.add(t);
            }
            Deque<Expr> st = new ArrayDeque<>();
            for (Token t : out){
                switch(t.type){
                    case NUM -> st.push(Expr.C(Double.parseDouble(t.text)));
                    case VAR -> st.push(new Var(t.text));
                    case FUNC -> { Expr a = st.pop(); st.push(applyFunc(t.text, a)); }
                    case OP -> { Expr b = st.pop(), a = st.pop(); st.push(applyOp(t.text, a, b)); }
                    default -> throw new RuntimeException("Bad token");
                }
            }
            if (st.size()!=1) throw new RuntimeException("Parse error");
            return st.pop();
        }
        private Token readNumber(){ int j=i; while(i<s.length() && (Character.isDigit(s.charAt(i))||s.charAt(i)=='.')) i++; return new Token(T.NUM, s.substring(j,i)); }
        private Token readIdentOrFunc(){ int j=i; while(i<s.length() && Character.isLetter(s.charAt(i))) i++; String id=s.substring(j,i);
            if (i<s.length() && s.charAt(i)=='(') return new Token(T.FUNC, id);
            return new Token(T.VAR, id);
        }
        private static int prec(String op){ return switch(op){ case "+","-"->1; case "*","/"->2; case "^"->3; default->0; }; }
        private static boolean rightAssoc(String op){ return op.equals("^"); }
        private static Expr applyFunc(String f, Expr a){
            return switch(f){
                case "sin"->Expr.sin(a);
                case "cos"->Expr.cos(a);
                case "exp"->Expr.exp(a);
                case "ln" -> new Ln(a);
                default-> throw new RuntimeException("Unknown func "+f);
            };
        }
        private static Expr applyOp(String op, Expr a, Expr b){
            return switch(op){
                case "+"->Expr.add(a,b); case "-"->Expr.sub(a,b); case "*"->Expr.mul(a,b);
                case "/"->Expr.div(a,b); case "^"->Expr.pow(a,b);
                default->throw new RuntimeException("Unknown op "+op);
            };
        }
        enum T { NUM, VAR, FUNC, OP, LP }
        record Token(T type, String text){}
    }

    public static final class AnalyticFunction implements Function1D {
        private final Expr expr;
        private final Expr dexpr;
        public AnalyticFunction(String expression, boolean enableSymbolic) {
            this.expr = new Parser(expression).parse();
            this.dexpr = enableSymbolic ? this.expr.d() : null;
        }
        public double f(double x, EvalContext ctx){ return expr.eval(x, ctx); }
        public Function1D symbolicDerivative(){ return (dexpr==null)? null : (xx, c) -> dexpr.eval(xx, c); }
        public String pretty(){ return expr.pretty(); }
    }

    // ===================== TABULATED FUNCTIONS =====================
    public record XY(double x, double y) implements Comparable<XY>{ public int compareTo(XY o){ return Double.compare(this.x, o.x); } }
    public record XYPair(XY left, XY right){}
    public interface TabulatedStorage { double xmin(); double xmax(); XYPair bracket(double x); int size(); }

    public static final class TreeSetStorage implements TabulatedStorage {
        private final TreeSet<XY> set;
        public TreeSetStorage(Collection<XY> points){ this.set = new TreeSet<>(points); }
        public double xmin(){ return set.first().x(); }
        public double xmax(){ return set.last().x(); }
        public int size(){ return set.size(); }
        public XYPair bracket(double x){
            XY q = new XY(x,0);
            XY floor = set.floor(q); XY ceil=set.ceiling(q);
            if (floor==null) floor=set.first(); if (ceil==null) ceil=set.last();
            return new XYPair(floor, ceil);
        }
    }

    public static final class TreeMapStorage implements TabulatedStorage {
        private final TreeMap<Double,Double> map;
        public TreeMapStorage(Map<Double,Double> pts){ this.map = new TreeMap<>(pts); }
        public double xmin(){ return map.firstKey(); }
        public double xmax(){ return map.lastKey(); }
        public int size(){ return map.size(); }
        public XYPair bracket(double x){
            Map.Entry<Double,Double> le = map.floorEntry(x), ge = map.ceilingEntry(x);
            if (le==null) le=map.firstEntry(); if (ge==null) ge=map.lastEntry();
            return new XYPair(new XY(le.getKey(), le.getValue()), new XY(ge.getKey(), ge.getValue()));
        }
    }

    public static final class TabulatedFunction implements Function1D {
        private final TabulatedStorage storage;
        public TabulatedFunction(TabulatedStorage storage){ this.storage = storage; }
        public double f(double x, EvalContext ctx){
            double xmin=storage.xmin(), xmax=storage.xmax();
            if (x<=xmin) return storage.bracket(xmin).left().y();
            if (x>=xmax) return storage.bracket(xmax).right().y();
            XYPair b=storage.bracket(x); double x0=b.left().x(), y0=b.left().y(); double x1=b.right().x(), y1=b.right().y();
            if (x1==x0) return y0;
            double t=(x-x0)/(x1-x0); return y0+t*(y1-y0);
        }
        public double xmin(EvalContext ctx){ return storage.xmin(); }
        public double xmax(EvalContext ctx){ return storage.xmax(); }
        public String pretty(){ return "Tabulated["+storage.size()+" pts]"; }
    }

    public static final class CSV {
        public static List<XY> loadToList(Path path) throws IOException {
            List<XY> out=new ArrayList<>();
            for(String line:Files.readAllLines(path)){ line=line.trim();
                if(line.isEmpty()||line.startsWith("#")) continue;
                String[] p=splitCSV(line); if(p.length<2) continue;
                out.add(new XY(Double.parseDouble(p[0]), Double.parseDouble(p[1])));
            }
            out.sort(Comparator.naturalOrder()); return out;
        }
        public static Map<Double,Double> loadToMap(Path path) throws IOException {
            TreeMap<Double,Double> m=new TreeMap<>();
            for(String line:Files.readAllLines(path)){ line=line.trim();
                if(line.isEmpty()||line.startsWith("#")) continue;
                String[] p=splitCSV(line); if(p.length<2) continue;
                m.put(Double.parseDouble(p[0]), Double.parseDouble(p[1]));
            } return m;
        }
        private static String[] splitCSV(String line){
            if(line.contains(",")) return line.split(",");
            if(line.contains(";")) return line.split(";");
            return line.split("\\s+");
        }
        public static void saveSinCSV(Path path,double start,double end,double step)throws IOException{
            try(BufferedWriter w=Files.newBufferedWriter(path)){
                w.write("# x,y\n");
                for(double x=start;x<=end+1e-12;x+=step)
                    w.write(String.format(Locale.US,"%.8f,%.12f%n",x,Math.sin(x)));
            }
        }
    }

    // ===================== COMPUTATION & OUTPUT =====================
    public static final class Runner {
        private final DifferentiationScheme diff; private final double xStart,xEnd,xStep,h; private final EvalContext ctx;
        public Runner(DifferentiationScheme d,double xStart,double xEnd,double xStep,double h,EvalContext ctx){this.diff=d;this.xStart=xStart;this.xEnd=xEnd;this.xStep=xStep;this.h=h;this.ctx=ctx;}
        public void computeToFile(Function1D f, Path out)throws IOException{
            try(BufferedWriter w=Files.newBufferedWriter(out)){
                w.write("# "+f.pretty()+"\n# x\tf(x)\tf'(x)\n");
                for(double x=xStart;x<=xEnd+1e-12;x+=xStep){
                    double fx=f.f(x,ctx); double dfx=diff.df(f,x,h,ctx);
                    w.write(String.format(Locale.US,"%.8f\t%.12f\t%.12f%n",x,fx,dfx));
                }
            }
        }
    }

    // ===================== SIMPLE SWING PLOTTER =====================
    public static final class PlotFrame extends JFrame {
        public PlotFrame(String title,Function1D f,Function1D df,EvalContext ctx,double xs,double xe){
            super(title); setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
            setSize(900,600); setLocationRelativeTo(null); add(new PlotPanel(f,df,ctx,xs,xe));
        }
    }
    public static final class PlotPanel extends JPanel {
        private final Function1D f,df; private final EvalContext ctx; private final double xs,xe;
        public PlotPanel(Function1D f,Function1D df,EvalContext ctx,double xs,double xe){this.f=f;this.df=df;this.ctx=ctx;this.xs=xs;this.xe=xe;}
        protected void paintComponent(Graphics g){super.paintComponent(g); Graphics2D g2=(Graphics2D)g;
            g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING,RenderingHints.VALUE_ANTIALIAS_ON);
            int w=getWidth(),h=getHeight(); int pad=50;
            double[] ys=sample(f,xs,xe,800), yds=sample(df,xs,xe,800);
            double ymin=min(ys,yds),ymax=max(ys,yds); if(ymax==ymin){ymax=ymin+1;}
            g2.setColor(Color.GRAY); g2.drawRect(pad,pad,w-2*pad,h-2*pad);
            if(xs<=0&&xe>=0) g2.drawLine(xTo(w,pad,xs,xe,0),pad,xTo(w,pad,xs,xe,0),h-pad);
            if(ymin<=0&&ymax>=0) g2.drawLine(pad,yTo(h,pad,ymin,ymax,0),w-pad,yTo(h,pad,ymin,ymax,0));
            g2.setStroke(new BasicStroke(2f));
            g2.setColor(new Color(0,90,200)); drawCurve(g2,ys,xs,xe,w,h,pad,ymin,ymax);
            g2.setColor(new Color(200,80,0)); drawCurve(g2,yds,xs,xe,w,h,pad,ymin,ymax);
            g2.setColor(Color.BLACK); g2.drawString("f(x)",pad+10,pad+15); g2.drawString("f'(x)",pad+60,pad+15);
        }
        private static void drawCurve(Graphics2D g2,double[] ys,double xs,double xe,int w,int h,int pad,double ymin,double ymax){
            Path2D path=new Path2D.Double(); int n=ys.length;
            for(int i=0;i<n;i++){ double x=xs+(xe-xs)*i/(n-1);
                int X=xTo(w,pad,xs,xe,x), Y=yTo(h,pad,ymin,ymax,ys[i]);
                if(i==0) path.moveTo(X,Y); else path.lineTo(X,Y);}
            g2.draw(path);
        }
        private static int xTo(int w,int pad,double xs,double xe,double x){return pad+(int)Math.round((x-xs)/(xe-xs)*(w-2*pad));}
        private static int yTo(int h,int pad,double ymin,double ymax,double y){return h-pad-(int)Math.round((y-ymin)/(ymax-ymin)*(h-2*pad));}
        private double[] sample(Function1D f,double xs,double xe,int n){double[]a=new double[n];
            for(int i=0;i<n;i++){double x=xs+(xe-xs)*i/(n-1);a[i]=f.f(x,ctx);} return a;}
        private static double min(double[]a,double[]b){double m=Double.POSITIVE_INFINITY;for(double v:a)m=Math.min(m,v);for(double v:b)m=Math.min(m,v);return m;}
        private static double max(double[]a,double[]b){double m=Double.NEGATIVE_INFINITY;for(double v:a)m=Math.max(m,v);for(double v:b)m=Math.max(m,v);return m;}
    }

    // ===================== UTIL =====================
    private static String format(double v){ DecimalFormat df=new DecimalFormat("0.########"); return df.format(v); }

    // ===================== MAIN =====================
    public static void main(String[] args)throws Exception{
        Locale.setDefault(Locale.US);
        double X_START=1.5,X_END=6.5,X_STEP=0.05,H=1e-3;
        DifferentiationScheme diff=new FivePointCentral();
        Map<String,Double> params=new HashMap<>(); params.put("a",1.0); EvalContext ctx=new EvalContext(params);
        Runner runner=new Runner(diff,X_START,X_END,X_STEP,H,ctx);

        // f1
        Function1D f1=new AnalyticFunction("exp(-x^2)*sin(x)",true);
        runner.computeToFile(f1,Paths.get("out_f1.txt"));
        // f2 with a
        for(double a:new double[]{0.5,1.0,1.5}){
            ctx.params.put("a",a);
            Function1D f2=new AnalyticFunction("exp(-a*x^2)*sin(x)",true);
            runner.computeToFile(f2,Paths.get(String.format(Locale.US,"out_f2_a%.1f.txt",a)));
        }
        // Tabulated sin(x)
        Path csv=Paths.get("tab_sin.csv");
        if(!Files.exists(csv)) CSV.saveSinCSV(csv,X_START,X_END,X_STEP);
        List<XY> ptsList=CSV.loadToList(csv); Function1D fTabSet=new TabulatedFunction(new TreeSetStorage(ptsList));
        runner.computeToFile(fTabSet,Paths.get("out_f3_tabulated_treeset.txt"));
        Map<Double,Double> ptsMap=CSV.loadToMap(csv); Function1D fTabMap=new TabulatedFunction(new TreeMapStorage(ptsMap));
        runner.computeToFile(fTabMap,Paths.get("out_f3_tabulated_treemap.txt"));

        // GUI demo
        ctx.params.put("a",1.0);
        Function1D fShow=new AnalyticFunction("exp(-a*x^2)*sin(x)",true);
        Function1D dfShow=(fShow.symbolicDerivative()!=null)?fShow.symbolicDerivative():(x,c)->diff.df(fShow,x,H,c);
        SwingUtilities.invokeLater(()->new PlotFrame("f(x) & f'(x) — a="+ctx.p("a",1.0),fShow,dfShow,ctx,X_START,X_END).setVisible(true));

        System.out.println("Done. Files generated: out_f1.txt, out_f2_a*.txt, out_f3_tabulated_*.txt, tab_sin.csv");
    }
}