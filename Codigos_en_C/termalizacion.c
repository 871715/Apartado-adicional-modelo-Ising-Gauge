#include "funciones_dinamica.h"
#include "funciones_red.h"
#include "random.h"
#include <string.h>

extern int xp[];
extern int yp[];
extern int zp[];
extern int xm[];
extern int ym[];
extern int zm[];

void guardar_parametros_comparacion(int FLAG_INI, int N_sweps_entre_medidas, int N_medidas, const char* metodo, int indice) {
    char folder_param[256];
    snprintf(folder_param, sizeof(folder_param), "Resultados_simulacion/TERMALIZACION/%.2f/PARAMETROS", beta);
    
    char filename[256];
    
    // Crear el archivo con índice específico para cada método
    if (strcmp(metodo, "Wilson") == 0) {
        snprintf(filename, sizeof(filename), "%s/Wilson_%d.txt", folder_param, indice);
    } else {
        snprintf(filename, sizeof(filename), "%s/Onn_%d.txt", folder_param, indice);
    }

    // Crear el archivo
    FILE* fparam = fopen(filename, "w");
    if (!fparam) {
        printf("No se pudo crear el archivo de parámetros %s\n", filename);
        return;
    }

    // Escribir los parámetros
    fprintf(fparam, "L\t%d\n", L);
    fprintf(fparam, "J\t%f\n", J);
    fprintf(fparam, "beta\t%f\n", beta);
    if(FLAG_INI==1){
        fprintf(fparam, "estado_inicial\tcold\n");
    }else{
        fprintf(fparam, "estado_inicial\thot\n");
    }
    fprintf(fparam, "N_sweps_entre_medidas\t%d\n", N_sweps_entre_medidas);
    fprintf(fparam, "N_medidas\t%d\n", N_medidas);
    fprintf(fparam, "Metodo\t%s\n", metodo);
    fprintf(fparam, "Indice_configuracion\t%d\n", indice);
    fprintf(fparam, "Tamanio_loop\t10x10\n");
    fprintf(fparam, "Observables\t");
    if (strcmp(metodo, "Wilson") == 0) {
        fprintf(fparam, "tiempo,energia_plaqueta,w10,magnetizacion\n");
    } else {
        fprintf(fparam, "tiempo,energia_plaqueta,O10,magnetizacion\n");
    }

    fclose(fparam);
    printf("Archivo de parámetros creado: %s\n", filename);
}

int main(){
    
    int s[3*L*L*L], plaquetas[3*L*L*L];
    double probabilidades[5];
    int FLAG_INI;
    int N_sweps_entre_med = 1;
    int N_medidas = 5000;

    printf("=========================================\n");
    printf("COMPARACIÓN TERMALIZACIÓN: WILSON vs O'ₙ - 10×10\n");
    printf("beta = %.2f, L = %d, Pasos = %d\n", beta, L, N_medidas);
    printf("=========================================\n");

    // ✅ PRIMERO: Buscar el mayor índice existente para no sobreescribir
    char folder_inicial[256];
    snprintf(folder_inicial, sizeof(folder_inicial), "Resultados_simulacion/TERMALIZACION/%.2f/CONFIGURACION_INICIAL", beta);
    
    int indice_base = 0;
    char filename_test[256];
    FILE* ftest;
    
    while (1) {
        snprintf(filename_test, sizeof(filename_test), "%s/I_%d.txt", folder_inicial, indice_base);
        ftest = fopen(filename_test, "r");
        if (ftest) {
            fclose(ftest);
            indice_base++;
        } else {
            break;
        }
    }
    
    printf("Comenzando desde índice: %d\n", indice_base);

    for(int FLAG_INI = 0; FLAG_INI < 2; FLAG_INI++){
        
        printf("\n=== CONFIGURACIÓN INICIAL: %s ===\n", 
               FLAG_INI == 0 ? "ALEATORIA (HOT)" : "ORDENADA (COLD)");

        for(int j = 0; j < 10; j++){
            int indice_actual = indice_base + j + FLAG_INI * 10;
            printf("\n--- EJECUCIÓN %d/10 (Índice: %d) ---\n", j + 1, indice_actual);
            
            // Inicializar con semilla única para esta configuración
            inicializa_PR(12345 + indice_actual);
            
            vector_cociente_prob(probabilidades);
            inicializa_vectores_de_vecinos();
            
            // Crear configuración inicial y guardar archivo
            crea_configuracionInicial(FLAG_INI, s);
            dame_plaquetas(s, plaquetas);
            
            // COPIAR configuración para ambos métodos
            int s_wilson[3*L*L*L], plaquetas_wilson[3*L*L*L];
            int s_onn[3*L*L*L], plaquetas_onn[3*L*L*L];
            memcpy(s_wilson, s, sizeof(int) * 3 * L * L * L);
            memcpy(plaquetas_wilson, plaquetas, sizeof(int) * 3 * L * L * L);
            memcpy(s_onn, s, sizeof(int) * 3 * L * L * L);
            memcpy(plaquetas_onn, plaquetas, sizeof(int) * 3 * L * L * L);
            
            printf("🔷 EJECUTANDO MÉTODO WILSON 10×10 (Índice: %d)...\n", indice_actual);
            guardar_parametros_comparacion(FLAG_INI, N_sweps_entre_med, N_medidas, "Wilson", indice_actual);
            dinamica_metropolis_w10(N_sweps_entre_med, N_medidas, probabilidades, s_wilson, plaquetas_wilson, indice_actual);
            
            printf("🔶 EJECUTANDO MÉTODO O'ₙ 10×10 (Índice: %d)...\n", indice_actual + 20);
            guardar_parametros_comparacion(FLAG_INI, N_sweps_entre_med, N_medidas, "Onn", indice_actual);
            dinamica_metropolis_O10(N_sweps_entre_med, N_medidas, probabilidades, s_onn, plaquetas_onn, indice_actual + 20);
            
            printf("✅ COMPARACIÓN %d/10 COMPLETADA\n", j + 1);
        }
    }

    printf("\n=========================================\n");
    printf("ANÁLISIS DE VENTANAS PARA TERMALIZACIÓN\n");
    printf("=========================================\n");

    int k_ini = 0, k_final = 19;
    int N_ventana = 5;
    int N_salto = 3 * 35;

    crear_ventanas(k_ini, k_final, N_ventana, N_salto);
    crear_media_global(k_ini, k_final);

    printf("\n=========================================\n");
    printf ("TERMALIZACIÓN COMPLETADA\n");
    printf("=========================================\n");
    printf("Resultados guardados en:\n");
    printf("- Wilson 10×10: Resultados_simulacion/TERMALIZACION/%.2f/EVOLUCION/\n", beta);
    printf("- O'ₙ 10×10:    Resultados_simulacion/TERMALIZACION/%.2f/EVOLUCION/\n", beta);
    printf("- Parámetros:    Resultados_simulacion/TERMALIZACION/%.2f/PARAMETROS/\n", beta);
    printf("- Ventanas:      Resultados_simulacion/TERMALIZACION/%.2f/VENTANAS/\n", beta);
    printf("\nPara analizar la termalización:\n");
    printf("1. Revisar archivos en VENTANAS/EVOLUCION/\n");
    printf("2. Analizar COMPATIBILIDAD_0_19.txt\n");
    printf("3. Comparar tiempos de termalización Wilson vs O'ₙ\n");
    printf("=========================================\n");

    return 0;
}