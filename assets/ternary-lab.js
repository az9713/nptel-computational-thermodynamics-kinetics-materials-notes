/* Original offline composition geometry; no equilibrium prediction. */
(() => {
  const root=document.getElementById('ternary-lab');
  if(!root)return;
  const svg=root.querySelector('svg');
  const ns='http://www.w3.org/2000/svg';
  const phases=[[.8,.15,.05],[.1,.8,.1],[.15,.15,.7]];
  const colors=['#7dd3fc','#f5a97f','#a6da95'];
  const xy=x=>[100+600*(x[1]+x[2]/2),570-600*Math.sqrt(3)*x[2]/2];
  const add=(tag,attrs,text)=>{const e=document.createElementNS(ns,tag);for(const [k,v] of Object.entries(attrs))e.setAttribute(k,v);if(text)e.textContent=text;svg.append(e);return e;};
  add('polygon',{points:'100,570 700,570 400,50.384758',fill:'none',stroke:'#bbc7dc','stroke-width':2});
  [[100,600,'A'],[700,600,'B'],[400,30,'C']].forEach(([x,y,t])=>add('text',{x,y,fill:'#e6edf3','text-anchor':'middle','font-size':24},t));
  const points=phases.map(xy);
  const triangles=points.map((_,i)=>add('polygon',{fill:colors[i],'fill-opacity':.16,stroke:'none'}));
  add('polygon',{points:points.map(p=>p.join(',')).join(' '),fill:'none',stroke:'#e6edf3','stroke-width':2});
  points.forEach((p,i)=>{add('circle',{cx:p[0],cy:p[1],r:8,fill:colors[i]});add('text',{x:p[0]+10,y:p[1]-12,fill:colors[i],'font-size':22},['alpha','beta','gamma'][i]);});
  const dot=add('circle',{r:9,fill:'#fff',stroke:'#101827','stroke-width':2});
  const inputs=[...root.querySelectorAll('input')];
  const result=root.querySelector('[data-result]');
  const amountLabels=[...root.querySelectorAll('[data-amount]')];
  function update(){
    const weights=inputs.map(e=>Number(e.value));
    const sum=weights.reduce((a,b)=>a+b,0);
    inputs.forEach((e,i)=>e.setAttribute('aria-valuetext',weights[i]+' relative amount units'));
    if(!sum){dot.setAttribute('visibility','hidden');triangles.forEach(e=>e.setAttribute('visibility','hidden'));result.textContent='All relative amounts are zero. Choose a positive amount to define a mixture.';amountLabels.forEach(e=>e.textContent='undefined');return;}
    const fractions=weights.map(w=>w/sum);
    const composition=[0,1,2].map(i=>fractions.reduce((s,f,j)=>s+f*phases[j][i],0));
    const pos=xy(composition);dot.setAttribute('cx',pos[0]);dot.setAttribute('cy',pos[1]);dot.setAttribute('visibility','visible');
    triangles.forEach((e,i)=>{e.setAttribute('points',[pos,points[(i+1)%3],points[(i+2)%3]].map(p=>p.join(',')).join(' '));e.setAttribute('visibility','visible');});
    amountLabels.forEach((e,i)=>e.textContent=(100*fractions[i]).toFixed(2)+'%');
    result.textContent='Overall atomic composition: A '+(100*composition[0]).toFixed(3)+'%, B '+(100*composition[1]).toFixed(3)+'%, C '+(100*composition[2]).toFixed(3)+'%. Phase compositions remain fixed.';
  }
  inputs.forEach(e=>e.addEventListener('input',update));update();
})();
