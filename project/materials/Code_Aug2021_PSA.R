## Title: Propensity Score Analysis Using R
## Date: 04.14.2021
## Code by: YM Baek & Rinseo Park
## Seoul R Meetup: April 2021
## [성향점수분석 주요 7단계 중 3~6단계에 대한 간단한 실습]

library("tidyverse") #데이터관리 및 사전처리, 시각화를 위해 
library("MatchIt") #매칭기법을 위해 

# LaLonde 데이터셋 소개
# 직업훈련 프로그램의 임금인상 효과(ATT)를 추정할 것 
?lalonde
# treat=0,1 그룹간 차이가 어떠한지?
lalonde %>%
  group_by(treat) %>%
  slice_head(n=5)

# 1. treat 변수만 투입한 경우
lalonde %>% group_by(treat) %>% summarize(mean(re78))
lm(re78~treat, lalonde) %>% summary()

# 2. covariate 변수들을 추가한 경우
lm(re78~., lalonde) %>% summary()
# formula 설정
myformula = as.formula("re78 ~ age + educ + race + married + nodegree + re74 + re75")
(myformulaT = update(myformula,treat~.))
(myformulaY = update(myformula,.~treat+.))

# 3. 성향점수를 사용한 경우
glm_fit = glm(myformulaT, lalonde, family=binomial(link="logit"))

# 3-1) 성향점수가중 기법
# 추정하고자 하는 처치효과의 종류(ATT,ATC,ATE)에 따라 가중치 부여
myd = lalonde %>% 
  mutate(
    ps=glm_fit$fitted.values, #추정확률 추출
    iptw=ifelse(treat==1,1,1/(1-ps))  #IPTW = ATT 계산용 가중치
  )
# IPTW 부여후 회귀분석 
lm(myformulaY, myd, weights=iptw) %>% summary()

# 3-2) 성향점수매칭 기법: 그리디 매칭
set.seed(12345) #매칭결과를 정확하게 반복하기 위해
# 1단계: matchit() 함수 
m_out=matchit(myformulaT,               #원인처치 확률 추정식
              myd,                      #데이터 
              distance = "glm",         #성향점수 추정방식: GLM
              link = "linear.logit",    #성향점수 형태: 선형로짓
              method="nearest",         #그리디 매칭, 최인접사례 매칭
              replace=TRUE)             #반복추출 허용 
#summary(m_out) #매칭결과 점검 
#plot(summary(m_out)) #공변량 균형성 점검

# 2단계: match.data() 함수 
dm_out=match.data(m_out)
#dm_out

# 3단계: 모수통계기법 기반 ATT 계산
lm(myformulaY, dm_out, weights=weights) %>% summary()

# 성향점수의 형태?
predict(glm_fit, type="link") %>% head() #로짓
predict(glm_fit, type="response") %>% head() #확률
# 간단한 비교
tibble(logit=predict(glm_fit, type="link"),
       prob=predict(glm_fit, type="response")) %>%
  rownames_to_column() %>%
  pivot_longer(-rowname,names_to="type") %>%
  ggplot(aes(x=value,fill=type,color=type))+
  geom_histogram(alpha=.3,binwidth=.1)+
  geom_vline(xintercept=0,lty="dashed")+
  geom_vline(xintercept=1,lty="dashed")+
  facet_wrap(~type,ncol=1)+
  theme_bw()

# 비모수통계기법을 사용할 경우 아래와 같이
library("Zelig")

# 1) zelig() 함수: 효과추정치 추정모형 설정 
z_psm_att=zelig(myformulaY,            #효과추정치 추정식 
                data=dm_out,           #데이터 
                model="ls",            #OLS 모형을 사용한다는 의미
                weights="weights",     #가중치를 부여한다는 의미
                cite=FALSE)            #젤리그 인용방식에 대한 안내문 숨김 
z_psm_att

# 2) setx() 함수: 추정모형을 적용할 독립변수의 조건상정 
x_psm_att0=setx(z_psm_att,treat=0,data=dm_out) #통제집단가정 상황
x_psm_att1=setx(z_psm_att,treat=1,data=dm_out) #처치집단가정 상황

# 3단계: 1단계와 2단계를 근거로 기댓값(expected value, ev) 시뮬레이션
# 디폴트는 1000회, 만약 5000번의 시뮬레이션을 원하는 경우 num=5000 옵션추가 
s_psm_att0=sim(z_psm_att,x_psm_att0)  #통제집단의 경우 기댓값
s_psm_att1=sim(z_psm_att,x_psm_att1)  #처치집단의 경우 기댓값

#ATT의 값을 추정한 후 95%신뢰구간 계산 
ev_att0=get_qi(s_psm_att0,"ev")  #1000번의 통제집단 기댓값들
ev_att1=get_qi(s_psm_att1,"ev")  #1000번의 처치집단 기댓값들
PSM_ATT=ev_att1-ev_att0          #1000번 추정된 ATT
quantile(PSM_ATT,p=c(0.025,0.50,0.975)) %>% round(3) #ATT의 95%신뢰구간 

# 3-3) 성향점수층화 기법
set.seed(12345) #매칭결과를 정확하게 반복하기 위해
# 1단계: matchit() 함수 
m_out=matchit(myformulaT,               #원인처치 확률 추정식
              myd,                      #데이터 
              distance = "glm",         #성향점수 추정방식: GLM
              link = "linear.logit",    #성향점수 형태: 선형로짓
              method="subclass",        #성향점수층화
              subclass=6)               #디폴트값 K=6 
#summary(m_out) #매칭결과 점검 
#plot(summary(m_out)) #공변량 균형성 점검

# 2단계: match.data() 함수 
dm_out=match.data(m_out)
#dm_out

# 3단계: 모수통계기법 기반 ATT 계산
lm(myformulaY, dm_out, weights=weights) %>% summary()