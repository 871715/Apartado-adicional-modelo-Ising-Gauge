#include "funciones_red.h"

// Definición de los arrays globales
int xp[L*L*L];
int yp[L*L*L];
int zp[L*L*L];
int xm[L*L*L];
int ym[L*L*L];
int zm[L*L*L];


void inicializa_vectores_de_vecinos(void) {

    for(int i = 0; i < L*L*L; i++) {
        xp[i] = 1;
        yp[i] = L;
        xm[i] = -1;
        ym[i] = -L;
        zp[i] = L*L;
        zm[i] = -L*L;
    }


    for(int i=0;i<L*L;i++){
        xp[L-1+L*i]=-(L-1);
        xm[L*i]=L-1;
    }

    for(int j=0;j<L;j++){        
        for(int k=0;k<L;k++){
            yp[L*(L-1)+L*L*j+k]=-L*(L-1);
            ym[L*L*j+k]=L*(L-1);
        }
    }

    for(int m=0;m<L;m++){
        for(int n=0;n<L;n++){
                zp[L*m+n+L*L*(L-1)]=-L*L*(L-1);
                zm[L*m+n]=L*L*(L-1);
        }
    }
    
}

int plaqueta_xy(int Nodo, int *aristas){
    int x, y, z;
    int s_nx = aristas[3*Nodo];
    int s_ny = aristas[3*Nodo+1];
    int s_n_mas_x_y = aristas[3*(Nodo + xp[Nodo]) + 1];
    int s_n_mas_y_x = aristas[3*(Nodo + yp[Nodo])];
    return s_nx * s_n_mas_x_y * s_ny * s_n_mas_y_x;
}

int plaqueta_xz(int Nodo, int *aristas){
    int x, y, z;
    int s_nx = aristas[3*Nodo];
    int s_nz = aristas[3*Nodo+2];
    int s_n_mas_x_z = aristas[3*(Nodo + xp[Nodo]) + 2];
    int s_n_mas_z_x = aristas[3*(Nodo + zp[Nodo])];
    return s_nx * s_n_mas_x_z * s_nz * s_n_mas_z_x;
}

int plaqueta_yz(int Nodo, int *aristas){
    int x, y, z;
    int s_ny = aristas[3*Nodo+1];
    int s_nz = aristas[3*Nodo+2];
    int s_n_mas_y_z = aristas[3*(Nodo + yp[Nodo]) + 2];
    int s_n_mas_z_y = aristas[3*(Nodo + zp[Nodo]) + 1];
    return s_ny * s_n_mas_y_z * s_nz * s_n_mas_z_y;
}

void dame_plaquetas(int *aristas, int *plaquetas){
    int V=L*L*L;
    for(int i=0;i<V;i++){
        plaquetas[3*i]=plaqueta_xy(i,aristas);
        plaquetas[3*i+1]=plaqueta_xz(i,aristas);
        plaquetas[3*i+2]=plaqueta_yz(i,aristas);
    }
}

void coordenadas_nodo(int Nodo, int *x, int *y, int *z){
    *z = Nodo/(L*L);
    *y = (Nodo - (*z)*L*L)/L;
    *x = Nodo - (*z)*L*L - (*y)*L;
}

double energia_normalizada(int *plaquetas){
    int V=L*L*L;
    double suma=0.0;
    for(int i=0;i<V;i++){
        suma=suma+plaquetas[3*i]+plaquetas[3*i+1]+plaquetas[3*i+2];
    }
    return (-J*suma);
}

void desviacion_estandar(int n, double datos[], double *media, double *desviacion){
    double suma=0.0;

    for(int i=0;i<n;i++){
        suma+=datos[i];
    }
    *media=suma/n;
    suma=0.0;

    
    for(int i=0;i<n;i++){
        suma+=(*media-datos[i])*(*media-datos[i]);
    }
    
    *desviacion=sqrt(suma/n);
}   

int magnetizacion(int *aristas){
    int V=3*L*L*L;
    int suma=0;
    for(int i=0;i<V;i++){
        suma=suma+aristas[i];
    }
    return suma;
}

int vecino_n_xp(int Nodo, int n){
    for(int i=0;i<n;i++){
        Nodo=Nodo+xp[Nodo];
    }
    return Nodo;
}

int vecino_n_yp(int Nodo, int n){
    for(int i=0;i<n;i++){
        Nodo=Nodo+yp[Nodo];
    }
    return Nodo;
}

int vecino_n_zp(int Nodo, int n){
    for(int i=0;i<n;i++){
        Nodo=Nodo+yp[Nodo];
    }
    return Nodo;
}


int un_loop_z(int Nodo_inicial, int *arista, int n){
    int lado_sur=1;
    int lado_oeste=1;
    int Nodo_sur=Nodo_inicial;
    int Nodo_oeste=Nodo_inicial;
    int x_sur, y_sur, z_sur;
    int x_oes, y_oes, z_oes;
    for(int i=0;i<n;i++){
        lado_sur=lado_sur*arista[3*(Nodo_sur)];
        lado_oeste=lado_oeste*arista[3*(Nodo_oeste)+1];
        Nodo_sur=Nodo_sur+xp[Nodo_sur];
        Nodo_oeste=Nodo_oeste+yp[Nodo_oeste];
    }
    int lado_norte=1;
    int lado_este=1;
    int Nodo_norte=Nodo_oeste;
    int Nodo_este=Nodo_sur;

    for(int i=0;i<n;i++){
        lado_norte=lado_norte*arista[3*(Nodo_norte)];
        lado_este=lado_este*arista[3*(Nodo_este)+1];
        Nodo_norte=Nodo_norte+xp[Nodo_norte];
        Nodo_este=Nodo_este+yp[Nodo_este];
    }
    return lado_sur*lado_oeste*lado_norte*lado_este;
}

int un_loop_y(int Nodo_inicial, int *arista, int n){
    int lado_sur=1;
    int lado_oeste=1;
    int Nodo_sur=Nodo_inicial;
    int Nodo_oeste=Nodo_inicial;
    for(int i=0;i<n;i++){
        lado_sur=lado_sur*arista[3*(Nodo_sur)+2];
        lado_oeste=lado_oeste*arista[3*(Nodo_oeste)];
        Nodo_sur=Nodo_sur+zp[Nodo_sur];
        Nodo_oeste=Nodo_oeste+xp[Nodo_oeste];
    }
    int lado_norte=1;
    int lado_este=1;
    int Nodo_norte=Nodo_oeste;
    int Nodo_este=Nodo_sur;

    for(int i=0;i<n;i++){
        lado_norte=lado_norte*arista[3*(Nodo_norte)+2];
        lado_este=lado_este*arista[3*(Nodo_este)];
        Nodo_norte=Nodo_norte+zp[Nodo_norte];
        Nodo_este=Nodo_este+xp[Nodo_este];
    }
    return lado_sur*lado_oeste*lado_norte*lado_este;
}

int un_loop_x(int Nodo_inicial, int *arista, int n){
    int lado_sur=1;
    int lado_oeste=1;
    int Nodo_sur=Nodo_inicial;
    int Nodo_oeste=Nodo_inicial;
    for(int i=0;i<n;i++){
        lado_sur=lado_sur*arista[3*(Nodo_sur)+1];
        lado_oeste=lado_oeste*arista[3*(Nodo_oeste)+2];
        Nodo_sur=Nodo_sur+yp[Nodo_sur];
        Nodo_oeste=Nodo_oeste+zp[Nodo_oeste];
    }
    int lado_norte=1;
    int lado_este=1;
    int Nodo_norte=Nodo_oeste;
    int Nodo_este=Nodo_sur;

    for(int i=0;i<n;i++){
        lado_norte=lado_norte*arista[3*(Nodo_norte)+1];
        lado_este=lado_este*arista[3*(Nodo_este)+2];
        Nodo_norte=Nodo_norte+yp[Nodo_norte];
        Nodo_este=Nodo_este+zp[Nodo_este];
    }
    return lado_sur*lado_oeste*lado_norte*lado_este;
}

double prom_Wilson_loops(int n, int *aristas){

    int promedio=0;
    int Nodo_actual_z=0;
    int Nodo_actual_y=0;
    int Nodo_actual_x=0;

    for(int i=0;i<L-1;i++){
        promedio=promedio+un_loop_z(Nodo_actual_z,aristas,n)+un_loop_z(Nodo_actual_y,aristas,n)+un_loop_z(Nodo_actual_x,aristas,n);
        Nodo_actual_z=Nodo_actual_z+zp[Nodo_actual_z];
        Nodo_actual_x=Nodo_actual_x+xp[Nodo_actual_x];
        Nodo_actual_y=Nodo_actual_y+yp[Nodo_actual_y];
    }
    return promedio/(3*L);
}


void inicializa_nodos_wilson( int n, int m, int nodos_wilson[][m][2]){
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            nodos_wilson[i][j][0]=generador_plano(0,L*L*L);
            nodos_wilson[i][j][1]=generador_plano(0,3);
        }
    }
}

void dame_wilsons_nn(int *aristas, int *wilsons, int n){
    int V=L*L*L;
    for(int i=0;i<V;i++){
        wilsons[3*i]=un_loop_x(i,aristas,n);
        wilsons[3*i+1]=un_loop_y(i,aristas,n);
        wilsons[3*i+2]=un_loop_z(i,aristas,n);
    }
}



/*
void dame_O_nn(int *aristas, double *O, int n){
    int potencias_2[4*n],combinaciones,iteraciones[4*n],plaquetas[3*L*L*L];
    int indices[4*n];
    potencias_2[0]=1;
    iteraciones[0]=0;

    double prob_wilson=0,funcion_particion=0,suma=0;

    for(int i=1;i<4*n;i++){
        potencias_2[i]=potencias_2[i-1]*2;
        iteraciones[i]=0;
    }
    combinaciones=potencias_2[4*n-1]*2;

    for(int j=0;j<L*L*L;j++){
        for(int k=0;k<3;k++){//direcciones
            indices_loop_nn(j,indices,n,k);
            suma=0;
            funcion_particion=0;
            for(int i=1;i<=combinaciones;i++){
                for(int m=0;m<4*n;m++){
                    iteraciones[m]++;
                }
                for(int l=0;l<4*n;l++){
                    if(iteraciones[l]==potencias_2[l]){
                        aristas[indices[l]]*=-1;
                        iteraciones[l]=0;
                    }
                }
                dame_plaquetas(aristas,plaquetas);
                prob_wilson=exp(-beta*energia_normalizada(plaquetas));
                funcion_particion+=prob_wilson;
                if(k==0){
                    suma+=un_loop_x(j,aristas,n)*prob_wilson;
                }else if(k==1){
                    suma+=un_loop_y(j,aristas,n)*prob_wilson;
                }else if(k=2){
                    suma+=un_loop_z(j,aristas,n)*prob_wilson;
                }
            }
            O[3*j+k]=suma/funcion_particion;
        }
    }
}
*/
//A ver si está bien esta versión corregida
void dame_O_nn(int *aristas, double *O, int n) {
    int V = L * L * L;
    double probabilidades[5];
    int muestras_terma;
    int muestras_medida;
    if (n<=3){
        muestras_medida=300;
        muestras_terma=200;
    }
    else if (n<=6){
        muestras_medida=200;
        muestras_terma=150;
    }
    else{
        muestras_medida=150;
        muestras_terma=100;
    }
    vector_cociente_prob(probabilidades);
    
    printf("🔄 Calculando O'_%d con MC local (%d+%d muestras)...\n", n, muestras_terma, muestras_medida);
    
    #pragma omp parallel for
    for (int idx = 0; idx < 3 * V; idx++) {
        int j = idx / 3;  // nodo
        int k = idx % 3;  // dirección
        
        // Crear copia local para este loop
        int *aristas_local = malloc(3 * V * sizeof(int));
        int *plaquetas_local = malloc(3 * V * sizeof(int));
        memcpy(aristas_local, aristas, 3 * V * sizeof(int));
        
        // Identificar aristas relevantes para este loop
        int indices[4 * n];
        indices_loop_nn(j, indices, n, k);
        int num_aristas_loop = 4 * n;
        
        double suma_wilson = 0.0;
        int aceptadas = 0;
        
        // 1. TERMALIZACIÓN LOCAL
        for (int paso = 0; paso < muestras_terma; paso++) {
            for (int a = 0; a < num_aristas_loop; a++) {
                int arista_idx = indices[a];
                int posiciones[4];
                posicion_plaquetas(arista_idx, posiciones);
                int index_prob = indice_cociente_prob(plaquetas_local, posiciones);
                
                double r = fran();
                if (r < probabilidades[index_prob]) {
                    // Aplicar flip
                    for (int i = 0; i < 4; i++) {
                        plaquetas_local[posiciones[i]] = -plaquetas_local[posiciones[i]];
                    }
                    aristas_local[arista_idx] = -aristas_local[arista_idx];
                    aceptadas++;
                }
            }
        }
        
        // 2. MEDICIÓN
        for (int paso = 0; paso < muestras_medida; paso++) {
            // Un paso de Metropolis local
            for (int a = 0; a < num_aristas_loop; a++) {
                int arista_idx = indices[a];
                int posiciones[4];
                posicion_plaquetas(arista_idx, posiciones);
                int index_prob = indice_cociente_prob(plaquetas_local, posiciones);
                
                double r = fran();
                if (r < probabilidades[index_prob]) {
                    for (int i = 0; i < 4; i++) {
                        plaquetas_local[posiciones[i]] = -plaquetas_local[posiciones[i]];
                    }
                    aristas_local[arista_idx] = -aristas_local[arista_idx];
                }
            }
            
            // Medir Wilson loop
            double wilson_val;
            if (k == 0) wilson_val = (double)un_loop_x(j, aristas_local, n);
            else if (k == 1) wilson_val = (double)un_loop_y(j, aristas_local, n);
            else wilson_val = (double)un_loop_z(j, aristas_local, n);
            
            suma_wilson += wilson_val;
        }
        
        O[idx] = suma_wilson / muestras_medida;
        
        // Verificar resultado
        if (isnan(O[idx]) || isinf(O[idx])) {
            printf("⚠️  O[%d] inválido = %f, forzando a 0\n", idx, O[idx]);
            O[idx] = 0.0;
        }
        
        free(aristas_local);
        free(plaquetas_local);
    }
    
    printf("✅ O'_%d calculado exitosamente\n", n);
}
void indices_loop_nn(int Nodo_inicial, int *vector, int n, int direccion){
    int nodo=Nodo_inicial;
    if(direccion==0){
    for(int i=0;i<n;i++){
        vector[i]=3*nodo+1;
        nodo+=yp[nodo];
    }
    for(int i=0;i<n;i++){
        vector[i+n]=3*nodo+2;
        nodo+=zp[nodo];
    }
    for(int i=0;i<n;i++){
        nodo+=ym[nodo];
        vector[i+2*n]=3*nodo+1;
    }
    for(int i=0;i<n;i++){
        nodo+=zm[nodo];
        vector[i+3*n]=3*nodo+2;
    }
    }else if(direccion==1){
        for(int i=0;i<n;i++){
        vector[i]=3*nodo;
        nodo+=xp[nodo];
    }
    for(int i=0;i<n;i++){
        vector[i+n]=3*nodo+2;
        nodo+=zp[nodo];
    }
    for(int i=0;i<n;i++){
        nodo+=xm[nodo];
        vector[i+2*n]=3*nodo;
    }
    for(int i=0;i<n;i++){
        nodo+=zm[nodo];
        vector[i+3*n]=3*nodo+2;
    }
    }else if(direccion==2){
        for(int i=0;i<n;i++){
        vector[i]=3*nodo+1;
        nodo+=yp[nodo];
    }
    for(int i=0;i<n;i++){
        vector[i+n]=3*nodo;
        nodo+=xp[nodo];
    }
    for(int i=0;i<n;i++){
        nodo+=ym[nodo];
        vector[i+2*n]=3*nodo+1;
    }
    for(int i=0;i<n;i++){
        nodo+=xm[nodo];
        vector[i+3*n]=3*nodo;
    }
    }
}

void calculo_promedios_onn(int *s, int n, int m, int n_pasos, int n_pasos_entre_mediciones, int n_termalizacion, double probabilidades[5]) {
    int aceptadas = 0;
    int plaquetas[3*L*L*L];
    dame_plaquetas(s, plaquetas);

    // Termalización
    N_pasos_metropolis(n_termalizacion, s, plaquetas, probabilidades, &aceptadas);

    // Reservar memoria para resultados
    double **resultados_onn = malloc(n * sizeof(double *));
    double **promedios_onn = malloc(n * sizeof(double *));
    
    for (int i = 0; i < n; i++) {
        resultados_onn[i] = malloc(m * sizeof(double));
        promedios_onn[i] = malloc(3 * sizeof(double)); // media, desv, error
    }

    // Inicializar
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < m; j++) {
            resultados_onn[i][j] = 0.0;
        }
    }

    // Dinámica y medición
    int num_medidas = n_pasos / n_pasos_entre_mediciones;
    
    for(int medida = 0; medida < num_medidas; medida++) {
        double *O_current = malloc(3*L*L*L * sizeof(double));
        
        for(int tamaño = 0; tamaño < n; tamaño++) {
            int n_val = tamaño + 1; // tamaños de 1 a n
            dame_O_nn(s, O_current, n_val);
            
            // Promedio sobre toda la red para este tamaño
            double suma = 0.0;
            for(int idx = 0; idx < 3*L*L*L; idx++) {
                suma += O_current[idx];
            }
            resultados_onn[tamaño][medida] = suma / (3*L*L*L);
        }
        
        free(O_current);
        N_pasos_metropolis(n_pasos_entre_mediciones, s, plaquetas, probabilidades, &aceptadas);
    }

    // Calcular estadísticas
    for(int i = 0; i < n; i++) {
        double media = 0.0, varianza = 0.0;
        
        for(int j = 0; j < num_medidas; j++) {
            media += resultados_onn[i][j];
        }
        media /= num_medidas;
        
        for(int j = 0; j < num_medidas; j++) {
            double diff = resultados_onn[i][j] - media;
            varianza += diff * diff;
        }
        varianza /= num_medidas;
        double desviacion = sqrt(varianza);
        double error = desviacion / sqrt(num_medidas);
        
        promedios_onn[i][0] = media;
        promedios_onn[i][1] = desviacion;
        promedios_onn[i][2] = error;
        
        printf("O_nn (n=%d): media = %f, desv = %f, error = %f\n", 
               i+1, media, desviacion, error);
    }

    // Guardar resultados (similar a calculo_promedios_wilson)
    char filename[256];
    FILE* file = fopen("Resultados_simulacion/PROMEDIOS_ONN/resultados.txt", "w");
    if (file) {
        for(int i = 0; i < n; i++) {
            fprintf(file, "%d\t%f\t%f\t%f\n", i+1, 
                    promedios_onn[i][0], promedios_onn[i][1], promedios_onn[i][2]);
        }
        fclose(file);
    }

    // Liberar memoria
    for(int i = 0; i < n; i++) {
        free(resultados_onn[i]);
        free(promedios_onn[i]);
    }
    free(resultados_onn);
    free(promedios_onn);
}
